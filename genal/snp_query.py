import aiohttp
import asyncio
import numpy as np
import nest_asyncio
from tqdm.asyncio import tqdm_asyncio

# Using nest_asyncio to allow execution in notebooks
nest_asyncio.apply()

# Main function to start the event loop and run the asynchronous query
def async_query_gwas_catalog(snps, p_threshold=5e-8, return_p=False, return_study=False, 
                             max_associations=None, timeout=100):
    loop = asyncio.get_event_loop()
    results_global, errors, timeouts = loop.run_until_complete(query_gwas_catalog_coroutine(snps, p_threshold, return_p, return_study, max_associations, timeout))
    return results_global, errors, timeouts


# Function to query GWAS Catalog API for SNP associations
async def query_gwas_catalog_coroutine(snps, p_threshold=5e-8, return_p=False, return_study=False, 
                                       max_associations=None, timeout=100):
    """
    Query the GWAS Catalog API for SNP associations.
    
    Parameters:
        snps (list): List of SNPs to query.
        p_threshold (float): P-value threshold for filtering associations.
        return_p (bool): Whether to return the P-value of the association.
        return_study (bool): Whether to return the study ID of the association.
        max_associations (int): Maximum number of associations to return for each SNP.
        timeout (int): Timeout for each query in seconds.

    Returns:
        results_global (dict): Dictionary storing the SNP (keys) and results for each SNP: a list of single strings or tuples
        errors (list): List storing SNP for which the GWAS Catalog could not be queried
        timeouts (list): List storing SNP for which the timeout was reached
    """
    
    results_global = {}  # Dictionary storing the SNP (keys) and results for each SNP: a list of single strings or tuples
    errors = []  # List storing SNP for which the GWAS Catalog could not be queried
    timeouts = [] # List storing SNP for which the timeout was reached

    async def fetch(session, url, timeout=timeout): 
        try:
            async with asyncio.timeout(timeout):
                async with session.get(url) as response:
                    if response.status == 200:
                        return await response.json()
                    return None
        except asyncio.TimeoutError:
            return "TIMEOUT"

    async def process_snp(session, snp):
        #print(f"Processing SNP {snp}")
        
        results_snp = []  # List storing the results for each association found for this SNP
        
        base_url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{snp}/associations?projection=associationBySnp"
        base_data = await fetch(session, base_url, timeout=timeout)
        
        if base_data == "TIMEOUT":
            timeouts.append(snp)
        elif base_data:
            i = 0
            # Process each association found for this SNP
            for assoc in base_data.get('_embedded', {}).get('associations', []):

                # If there are already max_associations, stop the loop
                if max_associations and i >= max_associations:
                    break
                i += 1

                pvalue = assoc.get("pvalue", np.nan)
                # If the pvalue of the association does not pass the threshold, the association is not processed further nor reported 
                if pvalue < p_threshold:
                    trait = assoc.get("efoTraits", [])[0].get("trait", "")
                    
                    # If the return_study flag is active: query the page containing the GWAS Catalog study ID
                    if return_study:
                        study_url = assoc.get("_links", {}).get("study", {}).get("href", {})
                        study_data = await fetch(session, study_url, timeout=timeout)
                        study_id = "TIMEOUT" if study_data == "TIMEOUT" else study_data.get("accessionId", "") if study_data else "Not found"
                    else:
                        study_id = None
                        
                    # Return a tuple or a string depending on the return flags
                    if return_p and return_study:
                        result_assoc = (trait, "{:.4g}".format(pvalue), study_id)
                    elif return_p:
                        result_assoc = (trait, "{:.4g}".format(pvalue))
                    elif return_study:
                        result_assoc = (trait, study_id)
                    else:
                        result_assoc = trait
                    results_snp.append(result_assoc)
                    
                else:
                    continue
                
            # Clean the associations depending on the flag
            # If the P-value and Study ID are not returned, display each trait only once
            if not return_p and not return_study:
                results_snp = list(set(results_snp))
            # If the P-value must be returned, return each trait once with the lowest p-value
            elif return_p and not return_study:
                min_trait = {}
                for trait, pvalue in results_snp:
                    if trait not in min_trait or pvalue < min_trait[trait]:
                        min_trait[trait] = pvalue
                results_snp = [(trait, min_trait[trait]) for trait in min_trait]
                
            results_global[snp] = results_snp
        else:
            errors.append(snp)
    
    async with aiohttp.ClientSession() as session:
        tasks = [process_snp(session, snp) for snp in snps]
        await tqdm_asyncio.gather(*tasks)
    
    # Exclude timeouts from errors
    #errors = [error for error in errors if error not in timeouts]

    return results_global, errors, timeouts
