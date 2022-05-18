## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Silvana Schneider <schneider.sil@gmail.com>'

New submission

Found the following (possibly) invalid DOIs:
  DOI: 10.1002/bimj.201800391
    From: DESCRIPTION
    Status: Service Unavailable
    Message: 503

Service was unavailable at that moment. DOI: 10.1002/bimj.201800391 checks out when manually checking.

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Resubmission
This is a resubmission. In this version I have:

* Omitted the examples of the functions we don't want to export

* Replaced "\dontrun{}" with "\donttest{}" in the examples as they can be executed but take longer than 5 seconds
