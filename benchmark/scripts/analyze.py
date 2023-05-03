import os
import math

def get_adder_graph_costs(line, word_size, idx_offset, unsigned_inputs):
    overshifts = []
    line = line.replace("Adder graph: ", "")
    line = line.replace("{{", "")
    line = line.replace("}}", "")
    if line == "{}":
        return 0, 0 # empty adder graph
    nodes = line.split("},{")
    num_nodes = 0
    num_FAs_total = 0
    coefficients = [1]
    FAs_nodes = {}
    for node in nodes:
        elements = node.split(",")
        num = elements[1]
        num = num.replace("[", "")
        num = num.replace("]", "")
        num = int(num)
        if num == 0 or num == 1:
            continue
        coefficients.append(num)
        scaled_num = num
        if idx_offset > 0:
            post_add_shift = int(elements[3])
            scaled_num = num * (2**post_add_shift)
        original_shift_a = int(elements[5+idx_offset])
        original_shift_b = int(elements[8+idx_offset])
        a = elements[3+idx_offset]
        a = a.replace("[", "")
        a = a.replace("]", "")
        a = int(a)
        w_a = word_size + math.ceil(math.log2(abs(a)))
        additional_shift_a = 0
        while a % 2 == 0:
            additional_shift_a += 1
            a /= 2
        b = elements[6+idx_offset]
        b = b.replace("[", "")
        b = b.replace("]", "")
        b = int(b)
        w_b = word_size + math.ceil(math.log2(abs(b)))
        sub = -b in coefficients
        additional_shift_b = 0
        while b % 2 == 0:
            additional_shift_b += 1
            b /= 2
        b_original = b
        if sub:
            b_original = -b
        can_copy_sign = (num >= 0 and a >= 0) or (num >= 0 and b_original >= 0) or (num < 0 and a < 0) or (num < 0 and b_original < 0)
        num_shift_a = additional_shift_a + original_shift_a
        num_shift_b = additional_shift_b + original_shift_b
        num_shift = max(num_shift_a, num_shift_b)
        num_FAs = math.ceil(math.log2(abs(scaled_num))) + word_size
        w_num = num_FAs
        MSB_LUT_removed = False
        if not sub: # a + (b << x) OR (a << x) + b
            num_FAs -= num_shift
            if num_shift_a >= w_b or num_shift_b >= w_a:
                num_FAs = 0
                MSB_LUT_removed = True
        elif original_shift_a == 0: # a - (b << x)
            num_FAs -= num_shift_b
            if num_shift_b >= w_a:
                num_FAs = 0
                MSB_LUT_removed = True
        if unsigned_inputs and original_shift_a + w_a < w_num and original_shift_b + w_b < w_num:
            num_FAs -= 1
            MSB_LUT_removed = True
        if (not MSB_LUT_removed) and can_copy_sign and (not unsigned_inputs):
            num_FAs -= 1
        num_nodes += 1
        num_FAs_total += num_FAs
    return num_nodes, num_FAs_total

def analyze_file(filepath, word_size, unsigned_inputs, idx_offset):
    costs = []
    with open(filepath, "r") as f:
        for line in f.readlines():
            line = line.replace("\n", "").replace("\r", "")
            if "Adder graph:" not in line:
                continue
            num_adders, num_bits = get_adder_graph_costs(line, word_size, idx_offset, unsigned_inputs)
            costs.append((num_adders, num_bits))
    if len(costs) < 1:
        costs = [(-1, -1)] # no solution was found
    return costs

def write_csv(csv_filepath, result_container):
    with open(csv_filepath, "w") as f:
        # write table head
        f.write("benchType;instanceName;wordSize;unsignedInputs;idx;numAdders;numBits\n") # tikz is a bit picky when using "_" in file headings => just use camel case instead
        # write data
        for bench_type, dict_1 in result_container.items():
            for (shift, FA), dict_2 in dict_1.items():
                for instance_name, dict_3 in dict_2.items():
                    for word_size, dict_4 in dict_3.items():
                        for unsigned_inputs, cost_list in dict_4.items():
                            if "enumerate" in bench_type:
                                write_to_file = cost_list # for enumeration experiments write all solutions to result file
                            else:
                                write_to_file = [cost_list[-1]] # for all "normal" experiments only write the "best" solution to result file
                            for idx, (num_adders, num_bits) in enumerate(write_to_file):
                                f.write(f"{bench_type};{instance_name};{word_size};{unsigned_inputs};{idx};{num_adders};{num_bits}\n")
                                

def main():
    word_size_list = [8, 16, 24, 32] # analyze results for these word sizes
    unsigned_inputs_list = [True, False] # analyze results for unsigned and signed inputs separately
    result_base_dir = "benchmark/results" # base path where all results are located
    csv_filepath = f"{result_base_dir}/results.csv"
    result_container = {}
    for base_dir_entry in os.scandir(result_base_dir):
        if not base_dir_entry.is_dir():
            continue
        bench_type = base_dir_entry.name
        result_container[bench_type] = {}
        for specific_dir_entry in os.scandir(f"{result_base_dir}/{bench_type}"):
            if not specific_dir_entry.is_dir():
                continue
            setting_elements = specific_dir_entry.name.split("_")
            if len(setting_elements) != 4:
                continue
            shift = int(setting_elements[1])
            FA = int(setting_elements[3])
            result_container[bench_type][(shift, FA)] = {}
            print(f"analyzing '{bench_type}/shift_{shift}_FA_{FA}' now...")
            for result_file_entry in os.scandir(f"{result_base_dir}/{bench_type}/shift_{shift}_FA_{FA}"):
                if not result_file_entry.is_file():
                    continue
                result_filename = result_file_entry.name
                if not result_filename.endswith(".txt"):
                    continue
                instance_name = result_filename.replace(".txt", "")
                complete_filepath = f"{result_base_dir}/{bench_type}/shift_{shift}_FA_{FA}/{instance_name}.txt"
                result_container[bench_type][(shift, FA)][instance_name] = {}
                for word_size in word_size_list:
                    result_container[bench_type][(shift, FA)][instance_name][word_size] = {}
                    for unsigned_inputs in unsigned_inputs_list:
                        idx_offset = shift # idx offset = 1 if shift is enabled, and 0 else => adder graph string notation is different in that case
                        result_container[bench_type][(shift, FA)][instance_name][word_size][unsigned_inputs] = analyze_file(complete_filepath, word_size, unsigned_inputs, idx_offset)
    print(f"finished analyzing files, writing results to csv file '{csv_filepath}' now...")
    write_csv(csv_filepath, result_container)

if __name__=='__main__':
    main()
