destination="~/Documents/$(basename $PWD)"
echo "Sending to: $destination @ $1"
rsync -az --delete --info=progress2 -e ssh ./ $1:$destination