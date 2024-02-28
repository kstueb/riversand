__version__ = "1.2.2"

# 24.02.2024: Changes to Riversand.add_samples(), .add_samples_from_dataframe():
# Fuction .add_samples() now accepts pd.DataFrame input.
# When reading from file or importing a DataFrame, the function now fils empty
# columns with default data.
# Columns 'lon', 'long', 'longitude', ... are now all recognized as topographic
# data (changes to params.all_keys)

# 25.02.2024: Implemented restandardization to 07KNSTD (Be) or KNSTD (Al).
# The function get_textline() now restandardizes. plot_poly_fit() restandardizes
# and plots restandardized data. Also change to naming convention for plotting
# multi-catchments.
# no unit tests for get_textline() and poly_E_results()

# 27.02.2024: Changes to plots.py; variable params._out_path
# Added a function self.restandardize(); no unit test
