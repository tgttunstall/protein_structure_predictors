import pandas as pd

from psp.mcsm_ppi import build_mutation_line, format_results


def test_build_mutation_line():
    row = {"chain": "A", "wt_aa": "R", "residue_number": "282", "mut_aa": "W"}
    assert build_mutation_line(row) == "A R282W"


def test_format_results(tmp_path):
    html = """
    <html><body>
    <table>
        <tr><th>Index</th><th>Mutation</th><th>ddG</th></tr>
        <tr><td>1</td><td>R282W</td><td>-1.23</td></tr>
    </table>
    </body></html>
    """
    result_file = tmp_path / "mcsm_ppi_job.html"
    result_file.write_text(html)

    output_csv = tmp_path / "out.csv"
    format_results(tmp_path, output_csv)

    df = pd.read_csv(output_csv)
    assert "Mutation" in df.columns
    assert "Index" not in df.columns
    assert df.loc[0, "Mutation"] == "R282W"
