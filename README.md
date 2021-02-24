# EscapeProbability
## A Mathematica Package for Astronomical Radiative Transfer Calculations

EscapeProbability extends the capabilities of Mathematica to solve the radiative transfer problem (in the astronomical context) using the escape probability approximation. It includes functions to numerically (iteratively) solve the involved equations of radiative transfer and statistical equilibrium. The relevant properties of the medium are taken from LAMBDA, the Leiden Atomic and Molecular Database (http://www.strw.leidenuniv.nl/~moldata/). The package includes functions to display the properties of many species, like energy levels, allowed transitions, collision rates and so on. This package has been inspired by the software RADEX (summarized in Van der Tak, F.F.S., Black, J.H., Schöier, F.L., Jansen, D.J., van Dishoeck, E.F. 2007, A&A 468, 627) which is an open source one-dimensional non-LTE radiative transfer code, that uses the escape probability formulation assuming an isothermal and homogeneous medium without large-scale velocity fields. (http://www.strw.leidenuniv.nl/~moldata/radex.html)
Escape Probability

One approach to reduce the complexity of the radiative transfer problem is to use approximations in order to decouple the radiative transfer equation from the statistical balance. We stated above that a photon that is emitted at a position x1 can be absorbed at another position x2 and influence the level population there. Assuming an homogeneous, isothermal medium one can derive the probability of an absorption while passing a column of absorber. This assumes that the medium along the flight path of the photon is in the same state as at position x1. This is a strong assumption, but justified under certain conditions. For example, consider a outward expanding spherical shell. The atoms at outer radii have a velocity relative to the atoms at the center. This velocity difference produces a Doppler shift in frequency which makes an absorption in the expanding shell unlikely because the profile functions drops rapidly for shifted frequencies. This Doppler shift makes the expanding shell virtually transparent for photons emitted at the center and thus decouples the radiation from the matter in the shell. This particular approximation is called Sobolev approximation of LVG (large velocity gradient).
EscapeProbability - Package

Below we give an example for a typical radiative transfer calculation which can be used to interpret radio-astronomical data. We wish to calculate the emission of an interstellar cloud of carbon monoxide with a given size, density and temperature. Carbon monoxide is the second most abundant molecule in interstellar space, after molecular hydrogen. It is very stable and requires low temperatures and densities in order to be able to emit radio- and submm radiation. Hydrogen is much more abundant than carbon making H2 much more abundant than CO (in the order of 30000 H atoms per C atom) but it is much more difficult to observe, since it requires extreme conditions to emit radiation. For this reason, CO is used as H2 proxy.  

Let's calculate, what CO intensities would be expected from a molecular cloud with the following parameters:

 * gas density (in this case number density of H2, since CO collidides dominantly with H2): 1×104 cm-3
 * CO column density (column density is the product of density times length, i.e. counting all CO molecules along a line of sight): 1×1014 cm-2
 * kinetic gas temperature: 50 K (which is relatively warm for a molecular cloud)
 * FWHM (Full Width Half Maximum) Doppler line width of the emission lines (the particles are moving and hence subject to Doppler shifts in frequency): 5 km/s

The command to perform the calculations is EscapeProbabilityRun. The result is shown below. The output is given in tabular form.

<img style="float: right;" src="https://github.com/Markusroellig/EscapeProbability/blob/master/EP1.jpg" width="100%">
 

You can see, that we expect the 115 GHz line to have a line strength of 31 mK in the line center. The total, line-integrated flux will be 0.166 K km/s. With this information we could now estimate how long a given radio telescope would have to observe the cloud in order to achieve a sufficient signal to noise ratio to detect this line.

The package is still in development, but I hope that I removed the most obvious bugs. You are free to download the package for private use. If you wish to use results from EscapeProbability in any publication please contact me to clarify how to refer to my work. Same goes if you plan to use the package commercially.

Please contact me if you encounter any bugs or have ideas for improvement. 

### Installation instructions

Unpack the contents of the ZIP file in the $UserBaseDirectory or $BaseDirectory. This will create a directory EscapeProbability conatining all package files, including the documentation.
