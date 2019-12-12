for func in matmul1 matmul2 matmul_dag1 matmul_dag2;
do
	 gobjdump --disassemble=_${func} -j .text -w complex_ops.o > ${func}.s
	 echo ${func}:
	 ./opcnt.sh ${func}.s
	 echo
done
