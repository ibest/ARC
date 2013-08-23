echo "Splits performed:"
grep -c "Split" log.txt
echo "Assemblies finished:"
grep -c "Assembly finished in" log.txt
echo "Assemblies failed:"
grep -c "Assembly failed after" log.txt
echo "Assemblies killed:"
grep -c "Assembly killed after" log.txt 
echo "Targets completed:"
grep -c "did not incorporate any more reads" log.txt
echo "Killed by repeat detection:"
grep -c "repetitive" log.txt

