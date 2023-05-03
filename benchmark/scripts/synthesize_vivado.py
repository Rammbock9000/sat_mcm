################################################################################
## Compilation and Result Generation script for Vivado
## This script is modified from the scipt of FloPoCo from Florent de Dinechin
################################################################################

import os
import sys
import re
import string
import argparse

def report(text):
    print("vivado_runsyn: " + text)


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='This is an helper script that launches Xilinx Vivado synthesis and extracts resource consumption and critical path information')
    parser.add_argument('-i', '--implement', action='store_true', help='Go all the way to implementation (default stops after synthesis)')
    parser.add_argument('-v', '--vhdl', help='VHDL file names (default flopoco.vhdl) - pass multiple files as colon-separated list, with the toplevel being the first file')
    parser.add_argument('-e', '--entity', help='Entity name (default is last entity of the toplevel VHDL file)')
    parser.add_argument('-t', '--target', help='Target name (default is read from the toplevel VHDL file)')
    parser.add_argument('-f', '--frequency', help='Objective frequency (default is read from the VHDL file)')
    parser.add_argument('-j', '--jobs', help='Number of threads that vivado is allowed to use for synthesis and implementation')
    parser.add_argument('-s', '--standard', help='VHDL standard (can be "VHDL" or "VHDL 2008", default is "VHDL")')
    parser.add_argument('-d', '--disable_flatten_hierarchy', action='store_true', help='Disable flattening hierarchy to disable optimizations beyond entity borders (default is "True")')

    options=parser.parse_args()
    
    disable_flatten = options.disable_flatten_hierarchy
    if disable_flatten == None:
        disable_flatten = True

    if options.implement==True:
        synthesis_only=False
    else:
        synthesis_only=True

    if options.vhdl==None:
        filenames_str = "flopoco.vhdl"
    else:
        filenames_str = options.vhdl
    
    found_vhdl = False
    found_vhd = False
    filenames = filenames_str.split(":")
    for filename in filenames:
        if ".vhdl" in filename:
            found_vhdl = True
        elif ".vhd" in filename:
            found_vhd = True
    if found_vhdl and found_vhd:
        ending = "both"
    elif found_vhdl:
        ending = ".vhdl"
    elif found_vhd:
        ending = ".vhd"
    else:
        print("corrupt vhdl filenames")
        exit()
    if len(filenames) < 1:
        print("did not pass any vhdl filenames")
        exit()
    toplevel_filename = filenames[0]
    
    if options.jobs==None:
        print("using 1 thread for synth+impl; use argument --jobs (-j) to define the number of available threads")
        jobs = 1
    else:
        try:
            jobs = int(options.jobs)
        except:
            jobs = 1
            print("failed to convert arguments for number of jobs to integer; using default value jobs=1")
    
    if options.standard == None:
        print("did not set vhdl standard, using default (VHDL) set to VHDL 2008 to enable newer standard")
        standard = "VHDL"
    elif options.standard=="VHDL 2008" or "2008":
        standard = "VHDL 2008"
    else:
        standard = "VHDL"
        if options.standard!="VHDL":
            print("unrecognized vhdl standard (" + options.standard + "), using default (VHDL)")

    entity = options.entity
    if entity == None:
        raise Exception("need toplevel entity name")
    target = options.target
    if target == None:
        raise Exception("need target FPGA name")
    frequency = options.frequency
    if frequency == None:
        raise Exception("need frequency")

    if target.lower()=="kintex7":
        part="xc7k70tfbv484-3"
    elif target.lower()=="virtex7":
      part="xc7v2000t-2gflg1925"
    elif target.lower()=="zynq7000":
        part="xc7z020clg484-1"
    elif target.lower()=="ultrascale":
        part="xcvu190-flga2577-2-e"
    elif target.lower()=="ultrascale2":
        part="xcvu13p-fhga2104-2-e"
    else:
        raise Exception("Target " + target + " not supported")

    #remove the path & extension:
    filenames_pure = []
    for filename in filenames:
        pos=0
        while pos != -1:
            last_pos = pos
            pos = filename.find("/",pos+1)

        if last_pos != 0:
            filename_pure = filename[last_pos+1:] #path removed
        else:
            filename_pure = filename

        pos = filename_pure.find(".")
        filename_pure = filename_pure[:pos] #extension removed
        filenames_pure.append(filename_pure)
    toplevel_filename_pure = filenames_pure[0]
    
    workdir=os.getcwd()+"/"+toplevel_filename_pure+"_"+target+"_"+frequency
    
    if os.path.exists(workdir):
    	report("Skipping "+workdir+" (already exists)")
    	exit()

    report("Working dir: "+workdir)
    report("Synthesizing "+filename)

    report("   entity:       " +  entity)
    report("   target:       " +  target)
    report("   frequency:    " +  frequency + " MHz")


    os.mkdir(workdir)
    filenames_abs = []
    for filename in filenames:
        filename_abs = os.path.abspath(filename)
        report("   design file:    " +  filename_abs)
        filenames_abs.append(filename_abs)

    os.chdir(workdir)
    

    # Create the clock.xdc file. For this we use the previous frequency,
    # but since I am too lazy to parse the VHDL to find inouts and output names, I copy the set_input_delay etc from the FloPoCo generated file
    # This file was created by FloPoCo to be used by the vivado_runsyn utility. Sorry to clutter your tmp.
    xdc_file_name = workdir+"/clock.xdc" 
    report("writing "+xdc_file_name)
    clock_file = open(xdc_file_name,"w")
    period = 1000.0/float(frequency) # the frequency is in MHz and the period in ns
    clock_file.write("create_clock -name clk -period " + str(period) + " [get_ports clk]\n")
    clock_file.close()

    project_name = entity
    tcl_script_name = os.path.join(workdir, project_name+".tcl")

    report("writing " +  tcl_script_name)
    tclscriptfile = open( tcl_script_name,"w")
    tclscriptfile.write("# Synthesis of " + entity + "\n")
    tclscriptfile.write("# Generated by FloPoCo's vivado-runsyn.py utility\n")
    tclscriptfile.write("create_project " + project_name + " -part " + part + "\n")
    
    for filename_abs in filenames_abs:
        tclscriptfile.write("add_files -norecurse " + filename_abs + "\n")
    if ending == ".vhdl":
        tclscriptfile.write("set_property FILE_TYPE {" + standard + "} [get_files *.vhdl]\n")
    elif ending == ".vhd":
        tclscriptfile.write("set_property FILE_TYPE {" + standard + "} [get_files *.vhd]\n")
    else:
        tclscriptfile.write("set_property FILE_TYPE {" + standard + "} [get_files *.vhdl]\n")
        tclscriptfile.write("set_property FILE_TYPE {" + standard + "} [get_files *.vhd]\n")
    tclscriptfile.write("read_xdc " + xdc_file_name + "\n")
    tclscriptfile.write("update_compile_order -fileset sources_1\n")
    tclscriptfile.write("update_compile_order -fileset sim_1\n")
    
    if disable_flatten:
        tclscriptfile.write("set_property STEPS.SYNTH_DESIGN.ARGS.FLATTEN_HIERARCHY none [get_runs synth_1]\n")

    # Reporting files
    utilization_report_file = workdir + "/"  + project_name + "_utilization_"
    timing_report_file = workdir + "/" + project_name + "_timing_"
    power_report_file = workdir + "/" + project_name + "_power_"

    if synthesis_only:
        utilization_report_file +="synth.rpt"
        timing_report_file +="synth.rpt"
        power_report_file +="synth.rpt"
    else:
        utilization_report_file +="placed.rpt"
        timing_report_file +="synth.rpt" #should be synth!
        power_report_file +="synth.rpt" #should be synth!

    if synthesis_only:
        result_name = "synth_1"
        # The following command is great because it remove the OBUFs but where to get the area?
        tclscriptfile.write("synth_design -mode out_of_context -jobs " + str(jobs) + " \n")
    else:
        result_name = "impl_1"
        tclscriptfile.write("launch_runs " + result_name + " -jobs " + str(jobs) + "\n")
        tclscriptfile.write("wait_on_run " + result_name + "\n")
        tclscriptfile.write("open_run " + result_name + " -name " + result_name + "\n")
    tclscriptfile.write("report_utilization  -file "+  utilization_report_file +"\n")

    # Timing report
    tclscriptfile.write("report_timing -file " + timing_report_file + " \n")

    tclscriptfile.write("report_power -file " + power_report_file + " \n")
        
    tclscriptfile.close()

    vivado_command = ("vivado -mode batch -source " + tcl_script_name)
    report("Running Vivado: " + vivado_command)
    os.system(vivado_command)

    report("cat " + utilization_report_file)
    os.system("cat " + utilization_report_file)
    report("cat " + timing_report_file)
    os.system("cat " + timing_report_file)
