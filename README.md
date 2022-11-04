# Zeolite_CIF_Parser
Automatically download and parse all of the available .CIF files and auxiliary framework information from the [IZA Database of Zeolite Structures](https://america.iza-structure.org/IZA-SC/ftc_table.php).

Please cite the original crystallographic data source if you wish to use the data in your own projects:\
Ch. Baerlocher and L.B. McCusker\
Database of Zeolite Structures: http://www.iza-structure.org/databases/

An example of the BEA framework atom coordinates plotted with `matplotlib`:

https://user-images.githubusercontent.com/92121568/200044145-467b8613-9622-4c2d-87aa-646e6f5302c7.mp4


## Code Components
- [IZA_CIFs_scrape](Code/IZA_CIFs_scrape.py): downloads all available CIF files from the IZA-SC Database of Zeolite Structures. **Example:** [CIF file for the MFI framework](CIFs/MFI.cif). 
- [IZA_framework_info_scrape](Code/IZA_framework_info_scrape.py): extracts framework-specific information from the IZA-SC Database of Zeolite Structures. **Example:** [auxiliary info for the MFI framework](FrameworkInfo/MFI_data.csv). 
- [crystal_CIF_to_XYZ](Code/crystal_CIF_to_XYZ.py): converts CIF files to XYZ data by parsing and applying the included symmetry operations to the atomic coordinates in the primitive unit cell. **Example:** [generated XYZ file for the MFI framework](XYZs/MFI.csv). 
- [bond_search](Code/bond_search.py): searches for Si-O bonds based on distance between atoms. Currently configured only for pure SiO<sub>2</sub> frameworks. **Example:** [Si-O bonds found for the MFI framework](FrameworkBonds/MFI.csv). 
- [cycles](Code/cycles.py): implements an efficient modified breadth-first search algrorithm to find chordless cycles of a zeolite framework graph. Reduces the cycle set to a set of rings. **Example:** [coordinates of small ring vertices for the MFI framework](Rings/MFI_rings.csv). 

## Instructions:
All of the available frameworks posted on the IZA website (as of Nov 1, 2022) have already been processed and the resulting files are available in this repository.

To reprocess the database (in case of future updates), run the individual scripts in the order they are introduced in [Code Components](#code-components) to generate corresponding files. **Note:** you will need to have dependencies listed in the [requirements](requirements.txt) file installed. To quickly install, run the following in your terminal:
```
py -m pip install -r requirements.txt
```
