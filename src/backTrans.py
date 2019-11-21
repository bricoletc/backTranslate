from Bio import SeqIO
import argparse, logging
import os
from runner import fromProtein, fromDNA


def CLI_parser():
    parser = argparse.ArgumentParser(
        description="Backtranslate amino acid sequences with pairwise distance constraints"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-i",
        "--input_file",
        dest="input_file",
        help="Path to fasta file containing amino acid or DNA sequences",
        type=str,
    )
    group.add_argument(
        "-s",
        "--sequence",
        dest="sequence",
        help="amino acid sequence passed on command-line",
        type=str,
    )

    parser.add_argument(
        "-m",
        "--mode",
        dest="mode",
        choices=["dna", "protein"],
        help="mode to run the tool on: in protein mode, samples backtranslations, in dna mode, uses the dna sequences as samples directly.",
        type=str,
        required=True,
    )

    parser.add_argument(
        "-d",
        "--min_distance",
        dest="min_dist",
        help="Minimum distance (Hamming, as fraction of distinct nucleotides) between all returned sequences",
        type=float,
        default=0.3,
    )

    parser.add_argument(
        "-n",
        "--num_samples",
        dest="num_samples",
        help="Maximum number of backtranslated DNA samples to produce",
        type=int,
        default=1000,
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        dest="output",
        help="An existing directory where the output will go. This is a fasta of the sampled sequence(s) and a stats file",
        type=str,
    )

    parser.add_argument(
        "--output_prefix",
        dest="out_prefix",
        help="prefix for output files",
        type=str,
        default="results",
    )

    parser.add_argument(
        "--stats_header",
        dest="stats",
        action="store_true",
        help="Prints the header for the stats file",
    )

    parser.add_argument(
        "--forbidden",
        help="File path to DNA sequences that cannot appear in the sample; one sequence per line.",
        type=str,
    )

    return parser


def validate_args(args, parser):

    if args.input_file is not None and not os.path.exists(args.input_file):
        logging.error(f"file {args.input_file} does not exist")
        # parser.print_help()
        exit(1)

    if args.output is not None and not os.path.exists(args.output):
        logging.error(f"directory {args.output} does not exist")
        exit(1)

    if args.forbidden is not None and not os.path.exists(args.forbidden):
        logging.error(f"file {args.forbidden} does not exist")
        exit(1)

    if not 0 <= args.min_dist <= 1:
        logging.error(f"min distance {args.min_dist} is not between 0 and 1")
        exit(1)

    if args.mode == "dna":
        if args.input_file is None:
            logging.error(
                f"Must provide an input file, not a sequence, when running in DNA mode."
            )
            exit(1)


def run_fromProtein(args, output_path, forbidden):
    if args.sequence is not None:
        fromProtein(
            args.sequence,
            args.num_samples,
            args.min_dist,
            seq_ID="No_ID",
            output_path=output_path,
            forbidden=forbidden,
        )

    else:
        num_records = 0
        this_output_path = output_path
        for seq_record in SeqIO.parse(args.input_file, "fasta"):
            if args.output is not None:
                this_output_path = f"{output_path}_{seq_record.id}"
            num_records += 1
            fromProtein(
                seq_record.seq,
                args.num_samples,
                args.min_dist,
                seq_ID=seq_record.id,
                output_path=this_output_path,
                forbidden=forbidden,
            )

        if num_records == 0:
            logging.warning(f"no fasta records were parsed from {args.input_file}")


def run_fromDNA(args, output_path, forbidden):
    fromDNA(
        args.input_file,
        args.min_dist,
        seq_ID=None,
        output_path=output_path,
        forbidden=forbidden,
    )


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(levelname)s: %(message)s", level=logging.INFO
    )

    parser = CLI_parser()
    args = parser.parse_args()
    if args.stats:
        header = "\t".join(BackTranslate.current_stats_header)
        logging.info(f"{header}")
        exit(0)
    validate_args(args, parser)
    output_path = (
        os.path.join(args.output, args.out_prefix) if args.output is not None else None
    )
    forbidden = []
    if args.forbidden is not None:
        with open(args.forbidden) as in_file:
            forbidden = in_file.readlines()
            forbidden = [el.strip().upper() for el in forbidden]

    if args.mode == "protein":
        run_fromProtein(args, output_path, forbidden)
    elif args.mode == "dna":
        run_fromDNA(args, output_path, forbidden)
