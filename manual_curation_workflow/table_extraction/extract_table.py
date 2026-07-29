from pathlib import Path
from paddlex import create_pipeline


# -----------------------------------------------------------------------------
# User-editable settings
# -----------------------------------------------------------------------------
IMAGE_FILE = Path("data_table_image.png")
OUTPUT_FOLDER = Path("table_output")
EXCEL_OUTPUT = OUTPUT_FOLDER / "data_table_output.xlsx"
HTML_OUTPUT = OUTPUT_FOLDER / "data_table_output.html"


def main() -> None:
    """Extract a table image and save only Excel and HTML outputs."""
    if not IMAGE_FILE.exists():
        raise FileNotFoundError(
            f"Image not found: {IMAGE_FILE.resolve()}"
        )

    OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

    print("Loading table-recognition pipeline...")
    pipeline = create_pipeline(
        pipeline="table_recognition_v2"
    )

    print(f"Processing: {IMAGE_FILE.name}")

    results = pipeline.predict(
        input=str(IMAGE_FILE),
        use_doc_orientation_classify=False,
        use_doc_unwarping=False,
    )

    found_result = False

    for result_number, result in enumerate(results, start=1):
        if result_number > 1:
            raise RuntimeError(
                "More than one extraction result was returned. "
                "This demo expects one input image and one table result."
            )

        found_result = True

        # Save only the reconstructed Excel and HTML tables.
        result.save_to_xlsx(str(EXCEL_OUTPUT))
        result.save_to_html(str(HTML_OUTPUT))

    if not found_result:
        raise RuntimeError("The pipeline returned no extraction results.")

    print(
        "\nFinished.\n"
        f"Excel output: {EXCEL_OUTPUT.resolve()}\n"
        f"HTML output: {HTML_OUTPUT.resolve()}"
    )


if __name__ == "__main__":
    main()
