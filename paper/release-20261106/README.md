# Release 20261106

This directory contains the frozen English (68-page) and Chinese (57-page) manuscripts, their 16 TeX source files, the bilingual source archive, and the original release README and SHA-256 checksum manifest.

The source and PDF bytes match the 20261106 release. Build instructions are in the bilingual source archive: pdfLaTeX for `en/`, XeLaTeX with the Windows fonts selected by ctex for `zh/`.

## Numerical supplement

The complete numerical supplement is distributed separately through [Zenodo](https://doi.org/10.5281/zenodo.23755453). Its archive is 225,601,793 bytes and is not committed to Git. [Download the numerical archive](https://zenodo.org/records/23755453/files/Theorem%20%281%2B1.75%29%20on%20the%20Goldbach%20Conjecture%20-%20Numerical%20Supplement%20-%2020261106.zip?download=1), place it in this directory, and check all release artifacts against `SHA256SUMS.txt`.

```sh
sha256sum --check SHA256SUMS.txt
```

Extract the numerical archive and run the seven commands in its `README.txt`, with assertions enabled and generated outputs outside the supplement. Appendix I identifies the inputs retained by each check and the certificate correspondence.

The mathematical claims, constants and numerical data are unchanged by this repository update. The saved finite checks have the scope stated in Appendix I; they are not a completed formalization of the analytic proof.
