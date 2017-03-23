from argparse import ArgumentParser
import yaml
from subprocess import check_call

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--build", default="hg38")
    parser.add_argument("yaml", 
        help="Path to bcbio-nextgen sample configuration YAML file.")
    parser.add_argument("blastdb",
        help="Path to BLAST nucleotide database.")
    args = parser.parse_args()
    blastdb = args.blastdb
    base_cmd = ("python tcrBiter_betav2.py --build {args.build} --r1 {r1} "
                " --r2 {r2} --blastdb {blastdb}")
    with open(args.yaml) as in_handle:
        summary = yaml.load(in_handle)
    details =  summary["details"]
    files = [x["files"] for x in details]
    r1files = [x[0] for x in files]
    r2files = [x[1] for x in files]
    names = [x["description"] for x in details]
    for r1, r2, name in zip(r1files, r2files, names):
        print "Processing %s." % name
        print base_cmd.format(**locals())
        check_call(base_cmd.format(**locals()), shell=True)
