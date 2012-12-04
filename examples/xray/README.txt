Single wavelength X-ray measurement.

The sample used in this data set is a magnetic thin film on glass,
measured with Cu K-alpha X-rays to obtain the structural parameters.

-----
    glass | FePt seed | 55:45 FePt | 80:20 NiFe | Pt cap
          |   2 nm    |   15 nm    |   50 nm    |  3 nm
-----

For a complete description of the system see [#ODonovan02]_

This model allows you to explore various initial conditions and search
ranges.

References
----------

.. [#ODonovan02]_:
    KV O'Donovan, JA Borchers, CF Majkrzak, O Hellwig and EE Fullerton (2002).
    "Pinpointing Chiral Structures with Front-Back Polarized Neutron
    Reflectometry".  Physical Review Letters 88(6) 067201(4)
    http://dx.doi.org/10.1103/PhysRevLett.88.067201


glass = SLD(name="glass", rho= 15.086, irho= 0.503)
seed  = SLD(name="seed",  rho=110.404, irho=13.769)
FePt  = SLD(name="FePt",  rho= 93.842, irho=10.455)
NiFe  = SLD(name="NiFe",  rho= 63.121, irho= 2.675)
cap   = SLD(name='cap',   rho= 86.431, irho=13.769)

sample = (glass(0,7.460)
          | seed(22.9417, 8.817)
          | FePt(146.576, 8.604)
          | NiFe(508.784, 12.736)
          | cap(31.8477, 10.715)
          | air)

