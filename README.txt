Instruction to use N-bodies integrator.

Code: LabPrac.cpp
Starting: ./LP

First, you have opportunity use some macroses in Lab Prac.cpp.
  FLOAT_PRESICION, DOUBLE_PRESICION, LONG_DOUBLE_PRESICION,ANY_PRESICION - obviously, switch precision of calculations.
  For any precision uses mpfr library, also you could using a key -b  only with this macros.
  With MEMORY_FREEING_BY_CURRENT, data will recording in RAM. With INSTANT_OUTPUT_TO_FILE, will instant recording in file(very slow).
  A key --freeing_memory_current only with this.

	
Input data:
  file_initial.dat
    Initial data for Caushy problem.
    In each string first three value - coordinates, next - velocities. Order of initial data must match order masses of bodies in file_mass.dat
  file_mass.dat
    Contains bodies masses.
  file_size.dat
    Contains T(integration end time), N(number of steps), M(number of bodies), dim(problem dimension * 2 ), bits(number of bits in data type), and pairs of ordinal numbers of bodies(in order match order in file_mass.dat and file_initial.dat), for which will calculate elements of orbits.
  file_encounters.dat
    Contains pairs of ordinal numbers of bodies(in order match order in file_mass.dat and file_initial.dat), which could collid, and distances, equaling sum of bodies radiuses.

Output data:
  res.dat
    Time, coordinates, velocities.
  orbital_elements.dat
    Time, semimajor axis, eccentricity, inclination, periapsis distance.
  distance_between_pairs.dat
    Time, distances.


Default dimension of data:
  Sun mass, astronomical units, years.
  G=4*pi^2.

All keys take priority over data files.

Keys:
-b Number bits in mpf_class precision. A key -b must be included first, because otherwise precision of t, h and other will be default 64 bits.
-G Gravitational constant(if problem have another dimension).
-N Number steps(if use controlled step, then t/N will minimal step).
-h Size of step(if use controlled step, then h will minimal step).
-M Number bodies.
-t Integration end time.
-p Step's decrease condition is abs((E(t)-E(t+h))/E(t)) > presicion_energy, where E(t) is energy integral in point t. Default equal 0. By dint of a key can be changed value precision_energy.
-d Problem dimension * 2.
--integrator Integrator name("RK4","RK5","LP").
--without_controlled_step 
--number_orbital_elements Number first pairs, for which do calculate orbital elements.
--number_encounters Number first pairs, for which do encounter test.
--focal_body Focal body is body relative to which will calculate coordinates and velocities other bodies.
--recording_frequency Frequency to which data recording to file.
--freeing_memory_current Number of steps after which memory will free and recording to file. Can be equal "default"(1.0e-6).

