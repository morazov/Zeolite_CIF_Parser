import numpy as np
import pandas as pd
import requests
import re

# This script downloads all available CIF files from the IZA-SC Database
# of Zeolite Structures.

# Please cite the crystallographic data source if you wish to use the data:
# Ch. Baerlocher and L.B. McCusker
# Database of Zeolite Structures: http://www.iza-structure.org/databases/

# Get the list of framework names from IZA database website
url_str = 'https://america.iza-structure.org/IZA-SC/ftc_table.php'
page = requests.get(url_str)
extracted_text = pd.read_html(page.text)
names = extracted_text[1].to_numpy().flatten()
names = filter(lambda v: v == v, names)
names = [re.sub("\*|\-", "", framework) for framework in list(names)]
names.sort()
print(names)

# Get and save the .cif files for each framework
for name in names:
    url_str = f'https://america.iza-structure.org/IZA-SC/cif/{name}.cif'
    file_request = requests.get(url_str, allow_redirects=True)
    open(f'./CIFs/{name}.cif', 'wb').write(file_request.content)
