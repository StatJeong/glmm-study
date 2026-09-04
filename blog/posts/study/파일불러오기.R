# 1) 패키지 설치  
install.packages("tidycensus")
library(tidycensus)


# 2) API 키 등록 (한 번만 하면 됨)  
census_api_key("YOUR_CENSUS_API_KEY", install = TRUE)

# 이후 R 세션 재시작 또는 아래로 .Renviron 읽기  
readRenviron("~/.Renviron")

# 예: 2022년 ACS 5-year PUMS에서 일부 변수(person level) 가져오기
pums_data <- get_pums(
  variables = c("AGEP","SEX","RAC1P","HINCP","PUMA","SERIALNO"), 
  state = "CA",    # 원하는 주
  survey = "acs5", 
  year = 2022
)

head(pums_data)
