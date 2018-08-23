# gmx_makecyclictop
This script makes a cyclic peptide or cyclic nucleotide topology file for gromacs out of the (specially prepared) linear molecule topology file

###Detailed instructions for the current version of the script

####(1) First step, take your cyclic peptide/nucleotide pdb file and imagine it’s linear (for a moment).

Add using text editor additional 2 additional residues to N and C termini.

The additional residues should be as follows. The one added in the beginning should be taken from the end, and the one added
to the end should be taken from the beginning.

for example, your cyclic peptide is   MKSY

you have to make a pdb file of peptide **Y**MKSY**M**

I highlighted the additional residues. The idea is : M has to follow after Y in the beginning (in the final cyclic petide; imagine the sequence going in circles), and M has to follow after Y in the end.

How you make the PDB file: Just COPY the Y into the beginning, and M into the the end.

Example: 

You original PDB file:


ATOM      1  N   MET A   1      82.620  48.430 -43.970  1.00  0.00           N  
ATOM      2  CA  MET A   1      81.410  47.750 -44.520  1.00  0.00           C  
ATOM      3  CB  MET A   1      81.670  46.260 -44.760  1.00  0.00           C  
ATOM      4  CG  MET A   1      80.870  45.660 -45.990  1.00  0.00           C  
ATOM      5  SD  MET A   1      79.050  45.500 -45.570  1.00  0.00           S  
ATOM      6  CE  MET A   1      78.680  45.540 -47.360  1.00  0.00           C  
ATOM      7  C   MET A   1      80.260  48.000 -43.630  1.00  0.00           C  
ATOM      8  O   MET A   1      79.720  49.090 -43.630  1.00  0.00           O  
ATOM      9  N   LYS A   2      79.800  47.030 -42.770  1.00  0.00           N  
ATOM     10  CA  LYS A   2      78.600  47.090 -41.990  1.00  0.00           C  
ATOM     11  CB  LYS A   2      78.420  45.710 -41.280  1.00  0.00           C  
ATOM     12  CG  LYS A   2      79.630  45.060 -40.670  1.00  0.00           C  
ATOM     13  CD  LYS A   2      79.390  43.670 -40.080  1.00  0.00           C  
ATOM     14  CE  LYS A   2      80.680  42.900 -39.920  1.00  0.00           C  
ATOM     15  NZ  LYS A   2      81.570  43.630 -38.990  1.00  0.00           N  
ATOM     16  C   LYS A   2      78.580  48.190 -40.960  1.00  0.00           C  
ATOM     17  O   LYS A   2      79.660  48.520 -40.400  1.00  0.00           O  
ATOM     18  N   SER A   3      77.400  48.800 -40.790  1.00  0.00           N  
ATOM     19  CA  SER A   3      77.280  50.040 -40.030  1.00  0.00           C  
ATOM     20  CB  SER A   3      76.560  51.140 -40.790  1.00  0.00           C  
ATOM     21  OG  SER A   3      77.380  51.610 -41.840  1.00  0.00           O  
ATOM     22  C   SER A   3      76.610  49.650 -38.690  1.00  0.00           C  
ATOM     23  O   SER A   3      75.510  49.110 -38.660  1.00  0.00           O  
ATOM     24  N   TYR A   4      77.300  49.890 -37.600  1.00  0.00           N  
ATOM     25  CA  TYR A   4      76.860  49.700 -36.270  1.00  0.00           C  
ATOM     26  CB  TYR A   4      78.030  49.390 -35.290  1.00  0.00           C  
ATOM     27  CG  TYR A   4      78.720  48.100 -35.660  1.00  0.00           C  
ATOM     28  CD1 TYR A   4      79.780  48.090 -36.570  1.00  0.00           C  
ATOM     29  CE1 TYR A   4      80.340  46.910 -36.930  1.00  0.00           C  
ATOM     30  CZ  TYR A   4      79.960  45.710 -36.360  1.00  0.00           C  
ATOM     31  OH  TYR A   4      80.500  44.450 -36.590  1.00  0.00           O  
ATOM     32  CD2 TYR A   4      78.220  46.890 -35.110  1.00  0.00           C  
ATOM     33  CE2 TYR A   4      78.840  45.700 -35.460  1.00  0.00           C  
ATOM     34  C   TYR A   4      76.150  50.920 -35.700  1.00  0.00           C  
ATOM     35  O   TYR A   4      76.640  52.050 -35.770  1.00  0.00           O  


New File, after copying the residues.

**ATOM     24  N   TYR A   4      77.300  49.890 -37.600  1.00  0.00           N  
ATOM     25  CA  TYR A   4      76.860  49.700 -36.270  1.00  0.00           C  
ATOM     26  CB  TYR A   4      78.030  49.390 -35.290  1.00  0.00           C  
ATOM     27  CG  TYR A   4      78.720  48.100 -35.660  1.00  0.00           C  
ATOM     28  CD1 TYR A   4      79.780  48.090 -36.570  1.00  0.00           C  
ATOM     29  CE1 TYR A   4      80.340  46.910 -36.930  1.00  0.00           C  
ATOM     30  CZ  TYR A   4      79.960  45.710 -36.360  1.00  0.00           C  
ATOM     31  OH  TYR A   4      80.500  44.450 -36.590  1.00  0.00           O  
ATOM     32  CD2 TYR A   4      78.220  46.890 -35.110  1.00  0.00           C  
ATOM     33  CE2 TYR A   4      78.840  45.700 -35.460  1.00  0.00           C  
ATOM     34  C   TYR A   4      76.150  50.920 -35.700  1.00  0.00           C  
ATOM     35  O   TYR A   4      76.640  52.050 -35.770  1.00  0.00           O**  
ATOM      1  N   MET A   1      82.620  48.430 -43.970  1.00  0.00           N  
ATOM      2  CA  MET A   1      81.410  47.750 -44.520  1.00  0.00           C  
ATOM      3  CB  MET A   1      81.670  46.260 -44.760  1.00  0.00           C  
ATOM      4  CG  MET A   1      80.870  45.660 -45.990  1.00  0.00           C  
ATOM      5  SD  MET A   1      79.050  45.500 -45.570  1.00  0.00           S  
ATOM      6  CE  MET A   1      78.680  45.540 -47.360  1.00  0.00           C  
ATOM      7  C   MET A   1      80.260  48.000 -43.630  1.00  0.00           C  
ATOM      8  O   MET A   1      79.720  49.090 -43.630  1.00  0.00           O  
ATOM      9  N   LYS A   2      79.800  47.030 -42.770  1.00  0.00           N  
ATOM     10  CA  LYS A   2      78.600  47.090 -41.990  1.00  0.00           C  
ATOM     11  CB  LYS A   2      78.420  45.710 -41.280  1.00  0.00           C  
ATOM     12  CG  LYS A   2      79.630  45.060 -40.670  1.00  0.00           C  
ATOM     13  CD  LYS A   2      79.390  43.670 -40.080  1.00  0.00           C  
ATOM     14  CE  LYS A   2      80.680  42.900 -39.920  1.00  0.00           C  
ATOM     15  NZ  LYS A   2      81.570  43.630 -38.990  1.00  0.00           N  
ATOM     16  C   LYS A   2      78.580  48.190 -40.960  1.00  0.00           C  
ATOM     17  O   LYS A   2      79.660  48.520 -40.400  1.00  0.00           O  
ATOM     18  N   SER A   3      77.400  48.800 -40.790  1.00  0.00           N  
ATOM     19  CA  SER A   3      77.280  50.040 -40.030  1.00  0.00           C  
ATOM     20  CB  SER A   3      76.560  51.140 -40.790  1.00  0.00           C  
ATOM     21  OG  SER A   3      77.380  51.610 -41.840  1.00  0.00           O  
ATOM     22  C   SER A   3      76.610  49.650 -38.690  1.00  0.00           C  
ATOM     23  O   SER A   3      75.510  49.110 -38.660  1.00  0.00           O  
ATOM     24  N   TYR A   4      77.300  49.890 -37.600  1.00  0.00           N  
ATOM     25  CA  TYR A   4      76.860  49.700 -36.270  1.00  0.00           C  
ATOM     26  CB  TYR A   4      78.030  49.390 -35.290  1.00  0.00           C  
ATOM     27  CG  TYR A   4      78.720  48.100 -35.660  1.00  0.00           C  
ATOM     28  CD1 TYR A   4      79.780  48.090 -36.570  1.00  0.00           C  
ATOM     29  CE1 TYR A   4      80.340  46.910 -36.930  1.00  0.00           C  
ATOM     30  CZ  TYR A   4      79.960  45.710 -36.360  1.00  0.00           C  
ATOM     31  OH  TYR A   4      80.500  44.450 -36.590  1.00  0.00           O  
ATOM     32  CD2 TYR A   4      78.220  46.890 -35.110  1.00  0.00           C  
ATOM     33  CE2 TYR A   4      78.840  45.700 -35.460  1.00  0.00           C  
ATOM     34  C   TYR A   4      76.150  50.920 -35.700  1.00  0.00           C  
ATOM     35  O   TYR A   4      76.640  52.050 -35.770  1.00  0.00           O  
**ATOM      1  N   MET A   1      82.620  48.430 -43.970  1.00  0.00           N  
ATOM      2  CA  MET A   1      81.410  47.750 -44.520  1.00  0.00           C  
ATOM      3  CB  MET A   1      81.670  46.260 -44.760  1.00  0.00           C  
ATOM      4  CG  MET A   1      80.870  45.660 -45.990  1.00  0.00           C  
ATOM      5  SD  MET A   1      79.050  45.500 -45.570  1.00  0.00           S  
ATOM      6  CE  MET A   1      78.680  45.540 -47.360  1.00  0.00           C  
ATOM      7  C   MET A   1      80.260  48.000 -43.630  1.00  0.00           C  
ATOM      8  O   MET A   1      79.720  49.090 -43.630  1.00  0.00           O**

I highlighted the copied residues. Notice how atom ordering and residue numbering got messed up.

Next step is renumbering residues (by hand) - important.
Renumber first and last residues by hand so that they follow one after another; the TYR will be 0 instead of 4, and
the last MET will be 5 instead of 4.

**ATOM     24  N   TYR A   0      77.300  49.890 -37.600  1.00  0.00           N  
ATOM     25  CA  TYR A   0      76.860  49.700 -36.270  1.00  0.00           C  
ATOM     26  CB  TYR A   0      78.030  49.390 -35.290  1.00  0.00           C  
ATOM     27  CG  TYR A   0      78.720  48.100 -35.660  1.00  0.00           C  
ATOM     28  CD1 TYR A   0      79.780  48.090 -36.570  1.00  0.00           C  
ATOM     29  CE1 TYR A   0      80.340  46.910 -36.930  1.00  0.00           C  
ATOM     30  CZ  TYR A   0      79.960  45.710 -36.360  1.00  0.00           C  
ATOM     31  OH  TYR A   0      80.500  44.450 -36.590  1.00  0.00           O  
ATOM     32  CD2 TYR A   0      78.220  46.890 -35.110  1.00  0.00           C  
ATOM     33  CE2 TYR A   0      78.840  45.700 -35.460  1.00  0.00           C  
ATOM     34  C   TYR A   0      76.150  50.920 -35.700  1.00  0.00           C  
ATOM     35  O   TYR A   0      76.640  52.050 -35.770  1.00  0.00           O**  
ATOM      1  N   MET A   1      82.620  48.430 -43.970  1.00  0.00           N  
ATOM      2  CA  MET A   1      81.410  47.750 -44.520  1.00  0.00           C  
ATOM      3  CB  MET A   1      81.670  46.260 -44.760  1.00  0.00           C  
ATOM      4  CG  MET A   1      80.870  45.660 -45.990  1.00  0.00           C  
ATOM      5  SD  MET A   1      79.050  45.500 -45.570  1.00  0.00           S  
ATOM      6  CE  MET A   1      78.680  45.540 -47.360  1.00  0.00           C  
ATOM      7  C   MET A   1      80.260  48.000 -43.630  1.00  0.00           C  
ATOM      8  O   MET A   1      79.720  49.090 -43.630  1.00  0.00           O  
ATOM      9  N   LYS A   2      79.800  47.030 -42.770  1.00  0.00           N  
ATOM     10  CA  LYS A   2      78.600  47.090 -41.990  1.00  0.00           C  
ATOM     11  CB  LYS A   2      78.420  45.710 -41.280  1.00  0.00           C  
ATOM     12  CG  LYS A   2      79.630  45.060 -40.670  1.00  0.00           C  
ATOM     13  CD  LYS A   2      79.390  43.670 -40.080  1.00  0.00           C  
ATOM     14  CE  LYS A   2      80.680  42.900 -39.920  1.00  0.00           C  
ATOM     15  NZ  LYS A   2      81.570  43.630 -38.990  1.00  0.00           N  
ATOM     16  C   LYS A   2      78.580  48.190 -40.960  1.00  0.00           C  
ATOM     17  O   LYS A   2      79.660  48.520 -40.400  1.00  0.00           O  
ATOM     18  N   SER A   3      77.400  48.800 -40.790  1.00  0.00           N  
ATOM     19  CA  SER A   3      77.280  50.040 -40.030  1.00  0.00           C  
ATOM     20  CB  SER A   3      76.560  51.140 -40.790  1.00  0.00           C  
ATOM     21  OG  SER A   3      77.380  51.610 -41.840  1.00  0.00           O  
ATOM     22  C   SER A   3      76.610  49.650 -38.690  1.00  0.00           C  
ATOM     23  O   SER A   3      75.510  49.110 -38.660  1.00  0.00           O  
ATOM     24  N   TYR A   4      77.300  49.890 -37.600  1.00  0.00           N  
ATOM     25  CA  TYR A   4      76.860  49.700 -36.270  1.00  0.00           C  
ATOM     26  CB  TYR A   4      78.030  49.390 -35.290  1.00  0.00           C  
ATOM     27  CG  TYR A   4      78.720  48.100 -35.660  1.00  0.00           C  
ATOM     28  CD1 TYR A   4      79.780  48.090 -36.570  1.00  0.00           C  
ATOM     29  CE1 TYR A   4      80.340  46.910 -36.930  1.00  0.00           C  
ATOM     30  CZ  TYR A   4      79.960  45.710 -36.360  1.00  0.00           C  
ATOM     31  OH  TYR A   4      80.500  44.450 -36.590  1.00  0.00           O  
ATOM     32  CD2 TYR A   4      78.220  46.890 -35.110  1.00  0.00           C  
ATOM     33  CE2 TYR A   4      78.840  45.700 -35.460  1.00  0.00           C  
ATOM     34  C   TYR A   4      76.150  50.920 -35.700  1.00  0.00           C  
ATOM     35  O   TYR A   4      76.640  52.050 -35.770  1.00  0.00           O  
**ATOM      1  N   MET A   5      82.620  48.430 -43.970  1.00  0.00           N  
ATOM      2  CA  MET A   5      81.410  47.750 -44.520  1.00  0.00           C  
ATOM      3  CB  MET A   5      81.670  46.260 -44.760  1.00  0.00           C  
ATOM      4  CG  MET A   5      80.870  45.660 -45.990  1.00  0.00           C  
ATOM      5  SD  MET A   5      79.050  45.500 -45.570  1.00  0.00           S  
ATOM      6  CE  MET A   5      78.680  45.540 -47.360  1.00  0.00           C  
ATOM      7  C   MET A   5      80.260  48.000 -43.630  1.00  0.00           C  
ATOM      8  O   MET A   5      79.720  49.090 -43.630  1.00  0.00           O**  

Notice how the residue numbers go from 0 to 5 now.

The atom numbering is in a mess but it’s not important, just the residue numbers are important.

####(2) next run your edited file (let's say *myfile.pdb*) through GROMACS using command something like:

*gmx pdb2gmx -f myfile.pdb -o myfile_gmx.pdb*

GROMACS creates topology file, write down it’s name (probably *topol.top*, but could be something else, like .itp file; if .itp is created, you should use it instead of .top for step 3).

####(3) Now run the perl script on that top file:

*gmx_makecyclictop.pl topfile*

where topfile is the topology file created in step (2).

The script creates new topology file, e.g. *topol_cyc.top*.

####(4) The hardest part, the topology file, is done. Now we just have to make sure we have the PDB (or GRO) file matching it.

For this, **delete the residues** added at step (1) **from the GROMACS generated file** (*myfile_gmx.pdb*)

Let's say the edited *myfile_gmx.pdb* without the first and last residues is named *myfile_gmx_edit.pdb*

Use that file for your simulation purposes, e.g. convert it to GRO and add simulation box:

*gmx editconf -f myfile_gmx_edit.pdb -o myfile_gmx_edit.gro -bt cubic -d 1.0*

Use *myfile_gmx_edit.gro* and *topol_cyc.top* for futher GROMACS simulations.

If you have further questions/bug reports, send email to visvaldas.kairys(at)bti.vu.lt.

Have fun!



