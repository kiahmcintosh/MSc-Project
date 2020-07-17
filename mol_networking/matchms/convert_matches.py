def convert_to_dictionary(scores):
    matches={}
    for score in scores:

        (reference, query, score, n_matching) = score
        if reference == query:
            continue
        reference=reference.metadata['scans']
        
        query=query.metadata['scans']

        if reference not in matches:
            matches[reference]={}
        matches[reference][query]={'cosine':score,'peaks':n_matching}
        if query not in matches:
            matches[query]={}
        matches[query][reference]={'cosine':score,'peaks':n_matching}
        
    
    return matches
