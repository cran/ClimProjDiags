# 0.3.0 (Release date: 2023-02-28)
- SelBox() and ShiftLon() to accept non-numerical data input  
- SelBox() uses the latitude and longitude dimension name instead of index  
- WeightedMean() uses multiApply::Apply inside  

# 0.2.1 (Release date: 2022-12-01)
- Fix the mistake that function "WeightedCells" was not included in the last submission.  

# 0.2.0 (Release date: 2022-11-04)
- New functions: ShiftLon, WeightedCells  
- Bugfix of Subset() when only one dimension left after subsetting and when 
parameter "indices" is not a list.
- Bugfix of WeightedMean() at the grid points that are across positive and
negative longitudes.


