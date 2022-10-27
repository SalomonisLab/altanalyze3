import pandas
import logging
import bioframe
from altanalyze3.utilities.helpers import TimeIt


def aggregate(args):
    with TimeIt():
        joined_df = None
        for int_location, jun_location, alias in zip(args.intcounts, args.juncounts, args.aliases):
            logging.info(f"Load counts from {int_location} and {jun_location} as {alias}")
            concatenated_df = pandas.concat(
                [
                    pandas.read_csv(
                        int_location,
                        usecols=[0, 1, 2, 4],
                        names=["chrom", "start", "end", alias],
                        index_col=["chrom", "start", "end"],
                        sep="\t",
                    ),
                    pandas.read_csv(
                        jun_location,
                        usecols=[0, 1, 2, 4],
                        names=["chrom", "start", "end", alias],
                        index_col=["chrom", "start", "end"],
                        sep="\t",
                    )
                ],
                axis=0
            )
            joined_df = concatenated_df if joined_df is None else joined_df.join(concatenated_df, how="outer")
        joined_df = joined_df.fillna(0)
        ref_df = pandas.read_csv(
            args.ref,
            usecols=[0, 1, 2, 3],
            names=["chrom", "start", "end", "name"],
            sep="\t",
        )
        closest_df = bioframe.closest(joined_df.index.to_frame(index=False), ref_df)
        print(closest_df)
        closest_df.set_index(["chrom", "start", "end"], inplace=True)
        closest_df.drop(labels=["chrom_", "start_", "end_", "distance"], axis=1, inplace=True)
        closest_df.rename(columns={"name_": "name"}, inplace=True)
        print(closest_df)
        joined_df = joined_df.join(closest_df, how="outer")
        print(joined_df)