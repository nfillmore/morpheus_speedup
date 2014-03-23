import argparse
import os
import subprocess
import sys

def get_morpheus_cl(build):
    morpheus_dir = os.path.dirname(os.path.abspath(__file__))
    morpheus_cl = os.path.join(morpheus_dir, "build", build, "morpheus_cl.exe")
    if not os.path.exists(morpheus_cl):
        print("The morpheus command, " + morpheus_cl + ", does not exist. Have "
              "you built morpheus yet?")
        sys.exit(1)
    return morpheus_cl

def run(morpheus_cl, morpheus_args, args):
    cmd = [args.mono, "--gc=" + args.gc]
    if args.profile is not None:
        cmd.append("--profile=" + args.profile)
    cmd.append(morpheus_cl)
    cmd.extend(morpheus_args)
    env = dict(os.environ)
    if args.gc_params is not None:
        gc_params = args.gc_params.replace("{nursery_size}",
                                           args.gc_nursery_size)
        env["MONO_GC_PARAMS"] = gc_params
    subprocess.check_call(cmd, env=env)

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="""
        This script runs the command-line version of Morpheus, morpheus_cl.exe,
        under Mono, using appropriate GC parameters by default. Command-line
        arguments other than those listed below are passed through to Morpheus.
    """, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--build", choices=("release", "debug"), default="release",
                   help="The build variant of morpheus_cl.exe to run.")
    p.add_argument("--mono", default="mono-sgen", help="The Mono executable.")
    p.add_argument("--gc", default="sgen", choices=("boehm", "sgen"),
                   help="The Mono garbage collector to use.")
    p.add_argument("--gc_params", default="major=marksweep-conc,minor=split,"
                                          "nursery-size={nursery_size}",
                   help="Parameters for Mono's garbage collector. The string "
                        "\"{nursery_size}\" will be replaced by the "
                        "--gc_nursery_size argument.")
    p.add_argument("--gc_nursery_size", default="32m",
                   help="The size of the Mono garbage collector's "
                        "\"nursery\", i.e., the pool for small, recent "
                        "objects that can be collected in parallel. This "
                        "should be set as large as possible without crashing.")
    p.add_argument("--profile",
                   help="Mono profiler options. If not given, don't profile.")
    args, unknown = p.parse_known_args()
    if len(unknown) == 0:
        p.print_help()
        sys.exit(1)
    morpheus_cl = get_morpheus_cl(args.build)
    run(morpheus_cl, unknown, args)
