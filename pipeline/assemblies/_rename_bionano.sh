#!/bin/sh

asm=$1
new_asm=${asm/.fasta/_rename.fasta}

mkdir -p map

asm_base=`basename $asm`
new_asm_base=`basename $new_asm`

if [ ! -e len/$asm_base.len ]
then
	echo "len/$asm_base.len not found."
	java -jar -Xmx2g /home/rhiea/codes/fastaContigSize.jar $asm
	mv $asm.len len/
fi

if [ -e map/$asm_base.list ]
then
	rm map/$asm_base.list
fi

if [ -e map/$asm_base.to.$new_asm_base.map ]
then
	rm map/$asm_base.to.$new_asm_base.map
fi

if [ -e map/$new_asm_base.list ]
then
	rm map/$new_asm_base.list
fi

if [ ! -e map/$asem_base.to.$new_asem_base.map ]
then
  for scaff in $(cut -f1 len/$asm_base.len);
  do
   scaff=`echo $scaff | awk '{print $1}'`;
   echo $scaff >> map/$asm_base.list;
   scaff2=${scaff/-/_};
   scaff2=${scaff2/:/_};
   echo -e $scaff"\t"$scaff2 >> map/$asm_base.to.$new_asm_base.map;
   echo $scaff2 >> map/$new_asm_base.list;
  done
fi
echo "java -jar -Xmx1g /home/rhiea/codes/fastaExtractFromList.jar $asm map/$asm_base.list $new_asm map/$new_asm_base.list"
java -jar -Xmx1g /home/rhiea/codes/fastaExtractFromList.jar $asm map/$asm_base.list $new_asm map/$new_asm_base.list

