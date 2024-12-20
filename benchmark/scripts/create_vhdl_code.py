import sys
import math
import os

def create_vhdl_code(mcm_str, input_word_size, output_filename):
    coeffs = []
    word_sizes = {1 : 0}
    indices = {1 : 0}
    output_shifts = {}
    inputs = {}
    signal_names = {1 : "x_1"}
    port_names = {}
    inst_names = {}
    subtract = {}
    copy_sign_left = {}
    copy_sign_right = {}
    nodes_str = mcm_str.split("},{")
    nodes_str_new = []
    # preprocess mcm string
    for i in range(len(nodes_str)):
        node = nodes_str[i]
        node = node.replace("{", "")
        node = node.replace("}", "")
        node = node.replace("'", "")
        node = node.replace("[", "")
        node = node.replace("]", "")
        node_elem = node.split(",")
        if len(node_elem) < 1:
            continue
        if node_elem[0] != "A":
            continue
        node = node.replace("A,", "")
        nodes_str_new.append(node)
    nodes_str = nodes_str_new
    # extract info from mcm string
    idx_counter = 0
    for node in nodes_str:
        elem = node.split(",")
        if len(elem) == 8:
            shift_o = 0
            idx_offset = 0
        elif len(elem) == 9:
            shift_o = int(elem[2])
            idx_offset = 1
        else:
            continue
        c = int(elem[0])
        if c < 0:
            signal_names[c] = f"x_m{abs(c)}"
            port_names[c] = f"y_m{abs(c)}"
            inst_names[c] = f"inst_m{abs(c)}_x"
        else:
            signal_names[c] = f"x_{c}"
            port_names[c] = f"y_{c}"
            inst_names[c] = f"inst_{c}_x"
        output_shifts[c] = shift_o
        coeffs.append(c)
        indices[c] = idx_counter
        idx_counter += 1
        word_sizes[c] = math.ceil(math.log2(abs(c)))
        c_in_left = int(elem[2+idx_offset])
        c_in_right = int(elem[5+idx_offset])
        if (c_in_right not in coeffs and -c_in_right in coeffs) or c_in_right == -1:
            subtract[c] = True
            c_in_right = -c_in_right
        else:
            subtract[c] = False
        shift_in_left = int(elem[4+idx_offset])
        shift_in_right = int(elem[7+idx_offset])
        inputs[c] = (signal_names[c_in_left], signal_names[c_in_right], shift_in_left, shift_in_right, word_sizes[c_in_left], word_sizes[c_in_right])
        if (c >= 0 and c_in_left >= 0) or (c < 0 and c_in_left < 0):
            copy_sign_left[c] = True
            copy_sign_right[c] = False
        elif (c >= 0 and c_in_right >= 0) or (c < 0 and c_in_right < 0):
            copy_sign_left[c] = False
            copy_sign_right[c] = True
        else:
            copy_sign_left[c] = False
            copy_sign_right[c] = False
        
    # create file from info
    with open(output_filename, "w") as f:
        f.write("\n")
        # libraries
        f.write("library ieee;\n")
        f.write("use ieee.std_logic_1164.all;\n")
        f.write("use ieee.numeric_std.all;\n")
        f.write("\n")
        # entity
        f.write("entity mcm is\n")
        f.write("  port (\n")
        for c in coeffs:
                f.write(f"    {port_names[c]} : out std_logic_vector({input_word_size}+{word_sizes[c]}-1 downto 0);\n")
        f.write(f"    x : in std_logic_vector({input_word_size}-1 downto 0)\n")
        f.write("  );\n")
        f.write("end mcm;\n")
        f.write("\n")
        # architecture head
        f.write("architecture mcm of mcm is\n")
        f.write(f"  signal x_1 : signed({input_word_size}-1 downto 0);\n")
        for c in coeffs:
            f.write(f"  signal {signal_names[c]} : signed({input_word_size}+{word_sizes[c]}-1 downto 0);\n")
        f.write("begin\n")
        f.write(f"  x_1 <= signed(x);\n")
        for c in coeffs:
            f.write(f"  {port_names[c]} <= std_logic_vector({signal_names[c]});\n")
        # instances
        f.write("\n")
        for c in coeffs:
            n_left = inputs[c][0]
            n_right = inputs[c][1]
            s_left = inputs[c][2]
            s_right = inputs[c][3]
            w_left = inputs[c][4]
            w_right = inputs[c][5]
            f.write(f"  {inst_names[c]} : entity work.add_node_v3\n")
            f.write("    generic map (\n")
            f.write(f"      w_x_i => {input_word_size + w_left},\n")
            f.write(f"      w_y_i => {input_word_size + w_right},\n")
            f.write(f"      w_o => {input_word_size + word_sizes[c]},\n")
            f.write(f"      s_x_i => {s_left},\n")
            f.write(f"      s_y_i => {s_right},\n")
            f.write(f"      s_o => {output_shifts[c]},\n")
            f.write(f"      copy_sign_x_i => {copy_sign_left[c]},\n".replace("T", "t").replace("F", "f"))
            f.write(f"      copy_sign_y_i => {copy_sign_right[c]},\n".replace("T", "t").replace("F", "f"))
            f.write(f"      sub => {subtract[c]}\n".replace("T", "t").replace("F", "f"))
            f.write("    )\n")
            f.write("    port map (\n")
            f.write(f"      x_i => {n_left},\n")
            f.write(f"      y_i => {n_right},\n")
            f.write(f"      z_o => {signal_names[c]}\n")
            f.write("    );\n")
        # architecture tail
        f.write("end mcm;\n")
    
    # create testbench from info
    if ".vhdl" in output_filename:
        testbench_filename = output_filename.replace(".vhdl", "_tb.vhdl")
    elif ".vhd" in output_filename:
        testbench_filename = output_filename.replace(".vhd", "_tb.vhd")
    else:
        # cannot create testbench because of weird naming scheme
        return
    with open(testbench_filename, "w") as f:
        # head
        f.write("\n")
        # libraries
        f.write("library ieee;\n")
        f.write("use ieee.std_logic_1164.all;\n")
        f.write("use ieee.numeric_std.all;\n")
        f.write("\n")
        # entity
        f.write("entity mcm_tb is\n")
        f.write("end mcm_tb;\n")
        f.write("\n")
        # architecture head
        f.write("architecture mcm_tb of mcm_tb is\n")
        f.write(f"  signal x_1 : std_logic_vector({input_word_size}-1 downto 0) := (others => '0');\n")
        f.write(f"  signal x_1_int : integer;\n")
        for c in coeffs:
            f.write(f"  signal {signal_names[c]} : std_logic_vector({input_word_size}+{word_sizes[c]}-1 downto 0);\n")
            f.write(f"  signal {signal_names[c]}_expected : integer;\n")
            f.write(f"  signal {signal_names[c]}_ok : std_logic;\n")
        f.write(f"  signal all_ok : std_logic;\n")
        f.write("begin\n")
        f.write(f"  all_ok <=")
        first = True
        for c in coeffs:
            if not first:
                f.write(" and")
            else:
                first = False
            f.write(f" {signal_names[c]}_ok")
        f.write(f";\n")
        f.write(f"  x_1_int <= to_integer(signed(x_1));\n")
        for c in coeffs:
            f.write(f"  {signal_names[c]}_expected <= {c} * x_1_int;\n")
            f.write(f"  {signal_names[c]}_ok <= '1' when {signal_names[c]}_expected = to_integer(signed({signal_names[c]})) else '0';\n")
        f.write(f"  DUT : entity work.mcm\n")
        f.write(f"    port map (\n")
        for c in coeffs:
            f.write(f"      {port_names[c]} => {signal_names[c]},\n")
        f.write(f"      x => x_1\n")
        f.write(f"    );\n")
        f.write(f"  process\n")
        f.write(f"  begin\n")
        f.write(f"    wait for 1ns;\n")
        f.write(f"    x_1 <= std_logic_vector(signed(x_1) + 1);\n")
        f.write(f"  end process;\n")
        f.write("end mcm_tb;\n")

def get_mcm_str(filepath):
    adder_graph = ""
    with open(filepath, "r") as f:
        for line in f.readlines():
            line = line.replace("\n", "")
            line = line.replace("\r", "")
            if "Adder graph:" not in line:
                continue
            line = line.replace("Adder", "")
            line = line.replace("graph", "")
            line = line.replace(":", "")
            line = line.replace(" ", "")
            adder_graph = line
    return adder_graph

def main():
    input_word_size = 8
    results_basepath = "benchmark/results"
    vhdl_basepath = "benchmark/vhdl"
    benchmarks = ["mcm", "unsigned_mcm"]
    sub_dir = "shift_0_FA_1"
    if not os.path.exists(f"{vhdl_basepath}/"):
        os.mkdir(f"{vhdl_basepath}/")
    for benchmark in benchmarks:
        if not os.path.exists(f"{vhdl_basepath}/{benchmark}/"):
            os.mkdir(f"{vhdl_basepath}/{benchmark}/")
        if not os.path.exists(f"{vhdl_basepath}/{benchmark}/{sub_dir}/"):
            os.mkdir(f"{vhdl_basepath}/{benchmark}/{sub_dir}/")
        complete_dir_path = f"{results_basepath}/{benchmark}/{sub_dir}/"
        for file_obj in os.listdir(os.fsencode(complete_dir_path)):
            filename = os.fsdecode(file_obj)
            mcm_str = get_mcm_str(f"{complete_dir_path}/{filename}")
            output_filename = f"{vhdl_basepath}/{benchmark}/{sub_dir}/{filename}".replace(".txt", ".vhd")
            create_vhdl_code(mcm_str, input_word_size, output_filename)

if __name__ == '__main__':
    main()
