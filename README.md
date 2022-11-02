# Zeolite_CIF_Parser
Automatically download and parse all of the available .CIF files and auxiliary framework information from the IZA Database of Zeolite Structures.

Please cite the original crystallographic data source if you wish to use the data in your own projects:\
Ch. Baerlocher and L.B. McCusker\
Database of Zeolite Structures: http://www.iza-structure.org/databases/

## Code Components
- [IZA_CIFs_scrape](): downloads all available CIF files from the IZA-SC Database of Zeolite Structures. **Example:** [CIF file for the MFI framework](). 
- [IZA_framework_info_scrape](): extracts framework-specific information from the IZA-SC Database of Zeolite Structures. **Example:** [auxiliary info for the MFI framework](). 
- [crystal_CIF_to_XYZ](): converts CIF files to XYZ data by parsing and applying the included symmetry operations to the atomic coordinates in the primitive unit cell. **Example:** [generated XYZ file for the MFI framework](). 
- [bond_search](): searches for Si-O bonds based on distance between atoms. Currently configured only for pure SiO<sub>2</sub> frameworks. **Example:** [Si-O bonds found for the MFI framework](). 
- [cycles](): implements an efficient modified breadth-first search algrorithm to find chordless cycles of a zeolite framework graph. Reduces the cycle set to a set of rings. **Example:** [coordinates of small ring vertices for the MFI framework](). 

## 
