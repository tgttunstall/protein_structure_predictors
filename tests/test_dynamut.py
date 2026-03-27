import pandas as pd

from psp.dynamut import build_mutation, format_results


def test_build_mutation():
    row = {"wt_aa": "F", "residue_number": "1174", "mut_aa": "C"}
    assert build_mutation(row) == "F1174C"


def test_format_results(tmp_path):
    html = """
    <html><body>
    <table>
        <tr><th>Mutation</th><th>ddG</th></tr>
        <tr><td>F1174C</td><td>0.12</td></tr>
    </table>
    </body></html>
    """
    result_file = tmp_path / "dynamut_job.html"
    result_file.write_text(html)

    output_csv = tmp_path / "out.csv"
    format_results(tmp_path, output_csv)

    df = pd.read_csv(output_csv)
    assert df.loc[0, "Mutation"] == "F1174C"
