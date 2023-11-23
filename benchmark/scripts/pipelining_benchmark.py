from start_benchmark import get_setup, do_it

def main():
    # tuple: (<bench-type>, <min-full-adders>, <allow-right-shift>, <pipelining>)
    # if <pipelining> is not specified (i.e. len(tuple) == 3), it is assumed to be 0
    bench_type_tuples = [("pmcm",0,1,1)]
    for bench_type_tuple in bench_type_tuples:
        do_it(bench_type_tuple, "kissat")
        print(f"finished benchmark {bench_type_tuple[0]} :)")

if __name__=='__main__':
    main()