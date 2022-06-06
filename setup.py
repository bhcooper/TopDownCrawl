import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name = 'TopDownCrawl',
    version = '1.0.1',
    author = 'Brendon Cooper',
    description = 'TopDownCrawl is a tool for aligning quantitative binding data for k-mers from experiments such as SELEX-seq and SMiLE-seq',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    url = "https://github.com/bhcooper/TopDownCrawl",
    packages=['TopDownCrawl'],
    install_requires = ['pandas', 'matplotlib', 'logomaker', 'fonts', 'font-roboto', 'xlrd', 'openpyxl'],
    entry_points = {
        "console_scripts": ['TopDownCrawl = TopDownCrawl.TopDownCrawl:main']
    },
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"],
)
