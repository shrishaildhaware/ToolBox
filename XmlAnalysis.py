
import xml.etree.ElementTree as ET

#-----------------------------------------------------------------
xmlStream = """<?xml version="1.0"?>
<data>
    <country name="Liechtenstein">
        <rank>1</rank>
        <year>2008</year>
        <gdppc>141100</gdppc>
        <neighbor name="Austria" direction="E"/>
        <neighbor name="Switzerland" direction="W"/>
    </country>
    <country name="Singapore">
        <rank>4</rank>
        <year>2011</year>
        <gdppc>59900</gdppc>
        <neighbor name="Malaysia" direction="N"/>
    </country>
    <country name="Panama">
        <rank>68</rank>
        <year>2011</year>
        <gdppc>13600</gdppc>
        <neighbor name="Costa Rica" direction="W"/>
        <neighbor name="Colombia" direction="E"/>
    </country>
</data>"""
#-----------------------------------------------------------------

root = ET.fromstring(xmlStream);
print(root.tag)

# Finding interesting elements
# Element has some useful methods that help iterate recursively over all the sub-tree below it (its children, their children, and so on). 
# For example, Element.iter():

for neighbor in root.iter('neighbor'):
    print(neighbor.attrib)
    
    
# Element.findall() finds only elements with a tag which are direct children of the current element.
# Element.find() finds the first child with a particular tag, and Element.text accesses the element’s text content. 
# Element.get() accesses the element’s attributes:

for country in root.findall('country'):
    rank = country.find('rank').text
    name = country.get('name')
    print(name, rank)
