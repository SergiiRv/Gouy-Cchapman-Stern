This is an implementation of the Gouy-Chapman-Stern (GChS)  model of potentials and ion distribution in nearest proximity and on neuronal membrane surface was done according to formulations provided by S. Genet & J. Burger [A] and M. Shapiro & F. Bezanilla [B]. The model was modified to take into account the temperature dependence of the ion’s hydrated radii according to S. Kuyucak [C]. Also, the model included the thermal dependence of the saline dielectric permittivity which was specifically measured for used in our experiments leech saline with the terahertz time-domain spectroscopy system and included concentration dependent component adapted and modified from A. Stogryn [D].

The model is implemented by S. Romanenko for the paper _________________________
Dr Sergii Romanenko; The University of Western Australia; CRAWLEY WA 6009;
+61 8 6488 7014; sergii.romanenko@uwa.edu.au; 

The theory and parametrs on the model is decribed in separate file SI3.docx

The model was inspired and substantially modified from:

A.	Genet, S., R. Costalat, and J. Burger, Acta Biotheoretica, 2000. 48(3/4): p. 273-287.
B.	Shapiro, M.G., et al., Infrared light excites cells by changing their electrical capacitance. Nature Communications, 2012. 3: p. 736.
C.	Kuyucak, S. and S.H. Chung, Temperature dependence of conductivity in electrolyte solutions and ionic channels of biological membranes. Biophys Chem, 1994. 52		(1): p. 15-24.
D.	Stogryn, A., Equations for Calculating the Dielectric Constant of Saline Water (Correspondence). IEEE Transactions on Microwave Theory and Techniques, 1971. 19		(8): p. 733-736.

Instructions:
1. Download files.
2. Unzip if needed. Place into the working folder.
3. Start the Matlab and redirect the path into the working library containing all the files. optimise the input parameters if needed.
4. run the GChS_run.m
----------------------------------
Additionally: The code saves the simulated data  in the excel file in the same folder.