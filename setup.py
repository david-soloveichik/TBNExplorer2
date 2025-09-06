from setuptools import find_packages, setup

setup(
    name="tbnexplorer2",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
    ],
    entry_points={
        "console_scripts": [
            "tbnexplorer2=tbnexplorer2.cli:main",
            "tbnexplorer2-filter=tbnexplorer2.filter_cli:main",
            "tbnexplorer2-ibot=extensions.ibot_cli:main",
        ],
    },
    include_package_data=True,
    data_files=[
        ('share/zsh/site-functions', [
            'completions/zsh/_tbnexplorer2',
            'completions/zsh/_tbnexplorer2-filter',
            'completions/zsh/_tbnexplorer2-ibot',
        ]),
    ],
    python_requires=">=3.8",
    author="TBN Explorer Team",
    description="A Python library and CLI tool for thermodynamics of binding networks (TBN) analysis",
    long_description="""A Python library and CLI tool for thermodynamics of binding networks (TBN) analysis""",
    long_description_content_type="text/markdown",
)
