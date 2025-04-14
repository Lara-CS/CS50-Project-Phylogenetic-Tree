import pytest
from main import hamming_formula, calculate_distance, read_seq
import numpy as np


def test_hamming_formula():
    assert hamming_formula("AAAAAAAAAAAA", "CAAAAAAAACCC") == 4
    with pytest.raises(ValueError):
        hamming_formula("AAAAAA", "AAAAAAAAAAAA")


def test_calculate_distance():
    a = calculate_distance(["AAAAAAAAAA", "AAAAAAAAAC", "AAAAAAAAGC", "AAAAAAACGC"])
    b = np.array([[0., 0.1, 0.2, 0.3], [0.1, 0., 0.1, 0.2], [0.2, 0.1, 0., 0.1], [0.3, 0.2, 0.1, 0.]])
    assert np.all(a == b)


def test_read_seq():
    assert str(read_seq("./test-data/test1.fna")) == 'AGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGC'
    with pytest.raises(TypeError):
        str(read_seq("test-data/test4.py"))


def test_main():
    with pytest.raises(FileNotFoundError):
        str(read_seq("./test-data/test5.fna"))


