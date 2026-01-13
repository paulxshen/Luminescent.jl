# Photonic PDK
Lumi Photonic PDK (process design kit) is a library of AI inverse designed photonic components for PIC (photonic integrated circuit). Unlike a traditional PDK, each component here is more akin to a reference design that shows what's possible. In most cases, they should be re-optimized for your specific needs and process node. **We can inverse design for any S-parameters objective with any geometry, wavelength, polarization, and layer stack.**

## Perfectly vertical grating coupler (PVGC)
Flexible, tolerant, perfectly vertical coupling between PIC and fiber / VCSEL / EEL / PD. Can inverse design for any modal profile. Sinble output design is compact and narrowband. Bidirectional design offers superior loss bandwidth tradeoff and is suited for broadband multiplexing applications. Below are sample designs for standard SMF-28 fiber mode at 1550nm wavelength.
### Single output PVGC
Loss: ~0.3dB (simulated)  
Bandwidth (1dB): ~10nm
### Narrowband Bidirectional PVGC
Loss: ~0.1dB (simulated)  
Bandwidth (1dB): ~10nm
### Broadband Bidirectional PVGC
Bandwidth (1dB): ~60nm
Bandwidth (3dB): ~100nm

![alt text](image.png)

## MMI aka splitters and combiners
MMIs with arbitrary split ratios and ports are easily inverse designed. 

## Wavelength domain multiplexer (WDM)
Compact wavelength domain multiplexer and demultiplexer with low insertion loss and crosstalk.

Updating... Please check back mid January