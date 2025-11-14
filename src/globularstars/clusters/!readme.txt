The folder "catalogues" contains results of astrometric classification of stars
in 170 globular clusters, using Gaia EDR3.
All stars in the given circular region with 5- or 6-parameter astrometry are included.
Columns are:
Gaia source_id;
celestial coordinates [deg]: ra, dec;
coordinates centered on the cluster [deg]: x,y;
parallax corrected for zero-point offset using the Lindegren et al. prescription [mas];
proper motion in ra, dec [mas/yr];
parallax uncertainty [mas];
PM uncertainty in ra, dec [mas/yr];
correlation coefficient between the two PM uncertainties;
G-band magnitude;
BP-RP;
source density [stars/arcmin^2], used to determine parallax/PM uncertainty scaling factors;
quality flag [bitfield - 4 possible combinations of two independent flags:
  lowest bit (0/1) denotes 6-parameter (0) or 5-parameter (1) solutions,
  next bit (0/2) distinguishes stars that passed all quality filters (2) from those that didn't (0)];
membership probability.

The folder "profiles" contains radial profiles of the rotational PM component and
the PM dispersion for these clusters.
Columns are:
distance from cluster centre [deg];
percentiles for the PM dispersion [mas/yr] (2.3, 15.9, 50, 84.1, 97.7 percentiles);
percentiles for the PM rotation   [mas/yr] (2.3, 15.9, 50, 84.1, 97.7 percentiles);
