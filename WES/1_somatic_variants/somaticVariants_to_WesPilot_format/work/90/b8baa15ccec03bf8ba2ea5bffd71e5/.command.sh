#!/usr/bin/env python
import pandas

filename = "intermediate_3_ready.xlsx"

processed_data = pandas.read_excel(filename)

# rename Chromosome, Position, End Position to Chr, Start and End
processed_data = processed_data.rename(columns={"Chromosome": "Chr",                                                 "Position": "Start",                                                 "End Position": "End"})
# write bed file
bed_file = processed_data[["Chr", "Start", "End"]]
bed_file.to_csv("output.bed", sep="	", index=False)
