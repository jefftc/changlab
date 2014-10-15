"""

Functions:
search_gene
find_antibodies


"""
# _search_gene_raw
# _parse_gene_search
# 
# _parse_gene_page_keyvalue
# _parse_gene_page_table
# _group_keyvalue_by_antibody

def search_gene(gene_name, wait=5):
    # Return a URL or None if gene not found.
    import urlparse

    html = _search_gene_raw(gene_name, wait=wait)
    for (name, url) in _parse_gene_search(html):
        if name.strip().upper() == gene_name.strip().upper():
            break
    else:
        return None
    
    x = urlparse.urljoin("http://www.antibodypedia.com", url)
    return x


def find_antibodies(gene_name, wait=5):
    # Return list of tuples (provider, antibody, clonality, western, elisa,
    # immunocytochemistry, immunoprecipitation, immunohistochemistry,
    # flow_cytometry).
    import urllib
    
    import timer
    
    # Search for the gene name.
    url = search_gene(gene_name, wait=wait)
    if not url:
        return []

    # Read and parse the antibodies.
    timer.wait(wait)
    html = urllib.urlopen(url).read()
    #open("test01.html", 'w').write(html)
    results_iter = _parse_gene_page_table(html)
    return list(results_iter)


def _search_gene_raw(gene_name, wait=5):
    # Searches for a gene and returns the raw HTML.  BUG: Will
    # generate Firefox browser in UI.
    from contextlib import closing
    #from selenium.webdriver import Firefox
    from selenium import webdriver
    from selenium.webdriver.support import expected_conditions as EC
    from selenium.webdriver.support.ui import WebDriverWait
    from selenium.webdriver.common.by import By

    import timer
    
    timer.wait(wait)
    #with closing(Firefox()) as browser:
    with closing(webdriver.PhantomJS()) as browser:
        URL = "http://www.antibodypedia.com/explore/%s" % gene_name
        browser.get(URL)

        # While searching, will show:
        # <div class="txtThree" id="search_results_title">
        #   Searching for '<strong>e2f1</strong>' ...</div>
        # 
        # When found, will show:
        # <div class="txtThree" id="search_results_title">
        #  Found <strong>8</strong> gene products for
        #  '<strong>e2f1</strong>'</div>
        wait = WebDriverWait(browser, 10)
        condition = EC.text_to_be_present_in_element(
            (By.ID, "search_results_title"), "Found")
        wait.until(condition)
        
        html = browser.page_source
        html = unicode(html).encode("utf-8")
    return html


def _parse_gene_search(html):
    # Yield (Name, URL)
    from bs4 import BeautifulSoup
    
    soup = BeautifulSoup(html)
    tables = soup.find_all("table", id="search_results_table")
    assert tables, "Could not find antibodies table"
    assert len(tables) == 1, "Found too many tables"
    results_table = tables[0]

    for row in results_table.find_all("tr"):
        if row.find_all("th"):
            # Skip the header or footer row.
            continue
        cols = row.find_all("td")
        # # Name Description Family Chromosome UniProt
        # Mouse_ortholog Antibodies
        if len(cols) == 1:
            col = cols[0]
            assert str(col.text).find("Search result is empty") >= 0
            continue
        assert len(cols) == 8, "Invalid cols: %s" % repr(cols)

        name = cols[1]
        hrefs = name.find_all("a")
        assert hrefs, "Could not find href"
        assert len(hrefs) == 1
        href = hrefs[0]
        name = str(href.text).strip()
        URL = href["href"]
        yield name, URL


def _parse_gene_page_table(html):
    # Yield tuples (provider, antibody, clonality, western, elisa,
    # immunocytochemistry, immunoprecipitation, immunohistochemistry,
    # flow_cytometry).

    for antibody_tups in _group_keyvalue_by_antibody(html):
        antibody_dict = {}
        for (key, value) in antibody_tups:
            antibody_dict[key] = value
        ad = antibody_dict
        x = ad["PROVIDER"], ad["ANTIBODY"], ad["CLONALITY"], \
            ad.get("WESTERN", ""), ad.get("ELISA", ""), \
            ad.get("IMMUNOCYTOCHEMISTRY", ""), \
            ad.get("IMMUNOPRECIPITATION", ""), \
            ad.get("IMMUNOHISTOCHEMISTRY", ""), \
            ad.get("FLOW_CYTOMETRY", "")
        yield x
    

def _parse_gene_page_keyvalue(html):
    # Yield tuples (KEY, VALUE).  Will have PROVIDER and
    # NUM_ANTIBODIES, followed by a number of ANTIBODY records.
    #
    # KEY                   VALUE
    # PROVIDER              <name>
    # NUM_ANTIBODIES        <num> (OPTIONAL)
    # ANTIBODY              <name>
    # CLONALITY             Polyclonal, Monoclonal
    # WESTERN               <evidence> (OPTIONAL)
    # ELISA                 <evidence> (OPTIONAL)
    # IMMUNOCYTOCHEMISTRY   <evidence> (OPTIONAL)
    # IMMUNOPRECIPITATION   <evidence> (OPTIONAL)
    # IMMUNOHISTOCHEMISTRY  <evidence> (OPTIONAL)
    # FLOW_CYTOMETRY        <evidence> (OPTIONAL)
    #
    # <EVIDENCE>
    # Supportive data in Antibodypedia
    # Data presented on provider website
    # Data in Antibodypedia (inconclusive)
    # Recommended by provider
    from bs4 import BeautifulSoup

    SCORE2DESC = {
        "/images/score_1.png" : "Supportive data in Antibodypedia",
        "/images/score_2.png" : "Data presented on provider website",
        "/images/score_3.png" : "Data in Antibodypedia (inconclusive)",
        "/images/score_4.png" : "Recommended by provider",
        }

    soup = BeautifulSoup(html)
    #print soup.prettify()
    tables = soup.find_all("table", id="antibodies_table")
    assert tables, "Could not find antibodies table"
    assert len(tables) == 1, "Found too many tables"
    antibodies_table = tables[0]
    #print antibodies_table.prettify()
    for row in antibodies_table.find_all("tr"):
        #print row.prettify(); continue
        if row.find_all("th"):
            # Skip the header or footer row.
            #print "HEADER OR FOOT"
            continue
        elif row.find_all("td", **{"class" : "provider"}):
            span = row.td.div.span.span
            descendants = list(span.descendants)
            assert len(descendants) == 6

            yield "PROVIDER", str(descendants[0].strip())
            
            # 26 antibodies
            x = descendants[4]
            x = x.replace("antibody", "")
            x = x.replace("antibodies", "")
            x = x.strip()
            x = int(x)  # make sure it's an int
            yield "NUM_ANTIBODIES", x
        elif row.find_all("td", colspan="10"):
            #<tr><td colspan="10" style="text-align:center;">
            #- No antibodies -
            #</td></tr>
            #
            # Must do this after "class=provider", because that also
            # has colspan=10.
            print row.prettify()
            assert str(row.text).find("No antibodies") >= 0
            continue
        elif row.find_all("td", **{"class" : "title"}):
            cols = row.find_all("td")
            assert len(cols) == 10
            #print row.prettify()

            yield "ANTIBODY", str(cols[1].text.strip())

            num_refs = 0
            x = cols[2].text.strip()
            if x:
                num_refs = int(x)
                yield "NUM_REFERENCES", num_refs
            
            x = str(cols[3].text.strip())
            # Multiple variants.  Don't bother checking.
            # Polyclonal
            # Monoclonal
            # <blank>
            # Polyclonal (Antigen purified)
            #assert x in ["Polyclonal", "Monoclonal", ""], x
            yield "CLONALITY", x

            c4, c5, c6, c7 = cols[4], cols[5], cols[6], cols[7]
            c8, c9 = cols[8], cols[9]

            if c4.find("img"):
                assert str(c4).find("Western blot") >= 0
                yield "WESTERN", SCORE2DESC[c4.img["src"]]
            
            if c5.find("img"):
                assert str(c5).find("ELISA") >= 0, str(cols[5])
                yield "ELISA", SCORE2DESC[c5.img["src"]]
            
            if c6.find("img"):
                assert str(c6).find("Immunocytochemistry") >= 0
                yield "IMMUNOCYTOCHEMISTRY", SCORE2DESC[c6.img["src"]]
                
            if c7.find("img"):
                assert str(c7).find("Immunoprecipitation") >= 0
                yield "IMMUNOPRECIPITATION", SCORE2DESC[c7.img["src"]]

            if c8.find("img"):
                assert str(c8).find("Immunohistochemistry") >= 0
                yield "IMMUNOHISTOCHEMISTRY", SCORE2DESC[c8.img["src"]]

            if c9.find("img"):
                assert str(c9).find("Flow cytometry") >= 0
                yield "FLOW_CYTOMETRY", SCORE2DESC[c9.img["src"]]
        else:
            raise AssertionError, "Unknown row\n%s" % row.prettify()


def _group_keyvalue_by_antibody(html):
    # Yield list of (key, value) associated with one antibody.
    provider = None
    antibody_tups = []  # list of (key, value)
    for x in _parse_gene_page_keyvalue(html):
        key, value = x
        # PROVIDER or ANTIBODY starts a new antibody.  Get rid of the
        # old information.
        if key in ["PROVIDER", "ANTIBODY"] and antibody_tups:
            # Don't store PROVIDER in the antibody_tups because there
            # are multiple ANTIBODY per PROVIDER.  Add it on now.
            assert provider
            antibody_tups = [("PROVIDER", provider)] + antibody_tups
            yield antibody_tups
            antibody_tups = []
        if key == "PROVIDER":
            provider = value
        elif key == "NUM_ANTIBODIES":
            pass
        else:
            antibody_tups.append((key, value))
    if antibody_tups:
        assert provider
        antibody_tups = [("PROVIDER", provider)] + antibody_tups
        yield antibody_tups
