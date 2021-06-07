import assaym


class args:
    pass
    

def test_deltag(capsys):
    args.assays="tests/cdc_assays.xlsx"
    args.sequences="tests/sequences.fasta"
    assaym.deltag(args)
    out,err=capsys.readouterr()
    with open('tests/expected.tsv', 'r') as file:
        expected = file.read()
    assert expected==out