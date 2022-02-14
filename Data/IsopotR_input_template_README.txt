save the excel spreadsheet as .csv

import in R via:
IsoplotR::read.data("file/to/your/data.csv", method = "...", format = x)

where method is one of 'U-Pb', 'Pb-Pb', 'Th-Pb', 'Ar-Ar', 'K-Ca', 'detritals', 'Rb-Sr', 'Sm-Nd', 'Re-Os', 'Th-U', 'U-Th-He', 'fissiontracks' or 'other'

and format is the number that is in the used excel-spreadsheet name (e.g. FT1 will be 1)

Look at ?IsoplotR::read.data if you need to know more