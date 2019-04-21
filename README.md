# modis

```python
sh preq.sh

python3 testgen.py 10 1000000 1000 > fasta/test_2.py

# in case you have sequences of different lengths
python3 msa.py .out fasta/test_2.py

python3 modis.py .out 10 10 20 0.3 1000 0.05
```
