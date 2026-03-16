echo "NODE annotaion" >  annotation; grep erpene *.anno | awk -F '[_.:]' '{print $5,$6,$7"."$8,$1,$2" TPS"}' OFS='_'  >> annotation; sed -i 's/ /\t/' annotation

