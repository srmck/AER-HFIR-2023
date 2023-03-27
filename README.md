# AER-HFIR-cycle-502

Simulations for the angle-encoded radiography experiment on the CG-4B beamline. Experiment scheduled for second half of cycle 502, which originally started March 28th 2023.
The cycle will be delayed for at least one month (as of (3/27/2023).

The sample will either be a 1mm borated aluminum diffraction grating or a Borated 3D printed mask with a slit size 2mm, 3mm, or 4mm.

### Issues fixed

1. Angle encoding is solved due to structure in empty beam.
2. Beam polarization can likewise be extracted from the empty beam run.
3. Image reconstruction can be done via the center guide field scan: the cosine and sine components are labeled by different currents.

### Top priority issues:

1. The image reconstruction algorithm needs to be automated, including the ptychographical phase stepper.
2. Ensure polarization transport after the last nutator to the s-bender spin analzyer.
3. Review the Anger camera data reduction script.

More details can be found at [AER Notes](https://docs.google.com/document/d/1wOUABEid8K96Qhme4a_Ta9n2TfiXdon4B7LQnB_AXkk/edit?usp=sharing).