# Cluster_stopping_calculation_drude_function

2020.10.28

## 概要

- Cluster_stopping_calculation_plasmon_pole_optical_limit と Cluster_stopping_calculation_plasmon_pole_with_dispersion に続き、三つ目に行った計算

- クラスター粒子を構成するイオンの電荷状態は、Kaneko(NIMB315(2013)76)を参考にした自己無撞着場計算で求めた(これは前の二つと同じ)
  - "0_Cx_average_charge.py"で前計算しておく
  - 結果は、resultsフォルダの"Cx_average_charge.json"に、単位原子当たりの入射エネルギー(keV)をキーとして保存するようにした

- Optical energy loss functionを、 Tan(Radiat.Environ.Biophys.45(2006)135)によって提案されているDrude type functionと仮定して計算を行った。すなわち Im{-1/ε(0,ω)} = <i>aω/((ω<sup>2</sup> - b<sup>2</sup>)<sup>2</sup>+(cω)<sup>2</sup>)</i>を用いた

- Energy loss functionは、Ashley(J.Phys.Condens.Matter3(1991)2741)のモデルを用いて、k > 0 の領域にOELFを拡張することで求めた

- ωの分散関係をちゃんと考えた計算。<i>ω<sub>k</sub> = ω<sub>p</sub> + (ℏk<sup>2</sup>/2m)</i>として計算を行った

- 波数についての積分は、<i>k</i> (not kappa!) を変数として行った
