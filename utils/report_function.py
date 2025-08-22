def report_function():
    html_content = """
    <html>
        <head><title>Report</title></head>
        <body>
            <h1>Sample Report</h1>
            <p>This is a generated report.</p>
        </body>
    </html>
    """

    with open("report.html", "w") as f:
        f.write(html_content)
    return "report.html"