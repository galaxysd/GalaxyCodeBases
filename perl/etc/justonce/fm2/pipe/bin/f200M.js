/** reads having a soft clip larger than 2 bases in beginning of read*/ 
function accept(rec) {   
  if (rec.getReadUnmappedFlag()) return false; 
  var cigar = rec.getCigar(); 
  if (cigar == null) return false;
  var readMatch = 0;
  for (var i=0;i < cigar.numCigarElements();++i) {
    var ce = cigar.getCigarElement(i);
    if (ce.getOperator().name() == "M") readMatch += ce.length;
  }
  if (readMatch > 200) return true;
}

accept(record); 
