# Map reads to various versions of KOLF assemblies

```sh
ver=v1.1

for pl in hifi ont_ul
do
  mkdir -p ${ver}.dip_$pl && cd ${ver}.dip_$pl
  ln -s ../${pl}.fofn input.fofn
  ln -s ../_submit.sh
  if [[ "$pl" == "hifi" ]]; then
    ./_submit.sh ../../assemblies/KOLF2.1J${ver}.fa ${ver}_dip.hifi map-pb
  elif [[ "$pl" == "ont_ul" ]]; then
    ./_submit.sh ../../assemblies/KOLF2.1J${ver}.fa ${ver}_dip.ont map-ont
  fi
  cd ../
done
```