import numpy as np
import pandas as pd
import requests
import re

# This script extracts framework-specific information from the IZA-SC Database
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

# Get the framework parameter data for each framework
for name in names:
    url_str = 'https://america.iza-structure.org/IZA-SC/'\
        f'framework.php?STC={name}'
    page = requests.get(url_str)
    extracted_text = pd.read_html(page.text)
    framework_info = pd.DataFrame(columns=['Parameter', 'Value'])
    framework_info.loc[0] = ["Crystal System", extracted_text[1].iloc[0, 2]]
    framework_info.loc[1] = ["Space Group", extracted_text[1].iloc[0, 3]]
    framework_info.loc[2] = ["Unit cell lengths",
                             ', '.join(str(e) for e in extracted_text[1]
                                       .iloc[1, 2:5].tolist())]
    framework_info.loc[3] = ["Unit cell angles",
                             ', '.join(str(e) for e in extracted_text[1]
                                       .iloc[2, 2:5].tolist())]
    framework_info.loc[4] = ["Volume", extracted_text[1].iloc[3, 3]]
    framework_info.loc[5] = ["Framework density", extracted_text[2].iloc[0, 2]]
    framework_info.loc[6] = ["Ring sizes (# T-atoms)",
                             extracted_text[2].iloc[2, 2]]
    framework_info.loc[7] = ["Channel dimensionality",
                             extracted_text[2].iloc[3, 2].split(':')[0]]
    framework_info.loc[8] = ["",
                             extracted_text[2].iloc[3, 2].split(':')[1]]
    framework_info.loc[9] = [
        "Maximum diameter of a sphere",
        ""]
    framework_info.loc[10] = [
        "that can be included",
        extracted_text[4].iloc[0, 2]]
    framework_info.loc[11] = [
        "Maximum diameter of a sphere",
        ""]
    framework_info.loc[12] = [
        "that can diffuse along",
        ', '.join(str(e) for e in extracted_text[5]
                  .iloc[0, 2:5].tolist())]

    framework_info.loc[13] = ["Accessible volume:",
                              extracted_text[6].iloc[0, 2]]
    framework_info['Value'] = framework_info['Value'].apply(
        lambda x: str(x).replace(u'\xa0', u''))
    framework_info['Value'] = framework_info['Value'].apply(
        lambda x: str(x).replace("Å3", "cubic Å"))
    framework_info['Value'] = framework_info['Value'].apply(
        lambda x: str(x).replace("Å", " Å"))
    framework_info['Value'] = framework_info['Value'].apply(
        lambda x: str(x).replace("  ", " "))
    framework_info.to_csv(
        f"./FrameworkInfo/{name}_data.csv", sep=';', index=False)
