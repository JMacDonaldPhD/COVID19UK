#' Download current COVID-19 Data (Public Health England)


currentNHSCaseData = function(){

  dat = read.csv(url("https://www.arcgis.com/sharing/rest/content/items/ca796627a2294c51926865748c4a56e8/data"),
                 stringsAsFactors = F)

  cases = as.numeric(dat$TotalCases)
  names(cases) = dat$NHSRNm
  cases
}

currentUACaseData = function(){
  dat = read.csv(url("https://www.arcgis.com/sharing/rest/content/items/b684319181f94875a6879bbc833ca3a6/data"),
           stringsAsFactors = F)

  cases = dat$TotalCases
  names(cases) = dat$GSS_NM
  cases
}

#' Work in Progress
pastUACaseData = function(date){

}
