from setuptools import setup, find_packages

setup(
    name="svjedi-tag",
    version="0.1.0",
    description="SVJedi-Tag: outil en ligne de commande pour analyse de variations structurales.",
    author="MTemperville",
    license="AGPL-3.0-only",
    packages=find_packages(),
    install_requires=[
        "networkx",
        "matplotlib",
        "biopython",
        "gfagraphs"
    ],
    entry_points={
        "console_scripts": [
            "svjedi-tag=svjedi_tag.svjedi-tag:main"
        ]
    },
    python_requires=">=3.11"
)

