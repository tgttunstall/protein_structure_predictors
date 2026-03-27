import pandas as pd

from psp.mcsm_lig import format_results


def test_format_results(tmp_path):
    html = """
    <html><body>
    <table>
        <tr><th>Mutation</th><th>Ligand</th></tr>
        <tr><td>D30N</td><td>065</td></tr>
    </table>
    </body></html>
    """
    result_file = tmp_path / "mcsm_lig_job.html"
    result_file.write_text(html)

    output_csv = tmp_path / "out.csv"
    format_results(tmp_path, output_csv)

    df = pd.read_csv(output_csv)
    assert df.loc[0, "Mutation"] == "D30N"
