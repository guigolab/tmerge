class MergedInfo():
    def __init__(self, hooks, speed = False, **kwargs):
        self.speed = speed

        hooks["transcript_added"].tap(self.add_attributes)
        hooks["transcripts_merged"].tap(self.add_dists_to_ends)
        hooks["transcripts_merged"].tap(self.update_contains)
        hooks["contig_merged"].tap(self.stringify)

    def stringify(self, transcripts):
        for transcript in transcripts:
            transcript.meta["3p_dists_to_3p"] = ",".join(transcript.meta["3p_dists_to_3p"])
            transcript.meta["5p_dists_to_5p"] = ",".join(transcript.meta["5p_dists_to_5p"])
            transcript.meta["contains"] = ",".join(transcript.meta["contains"])
            transcript.meta["contains_count"] = str(transcript.meta["contains_count"])

    def add_attributes(self, transcript):
        transcript.meta["3p_dists_to_3p"] = ["0"]
        transcript.meta["5p_dists_to_5p"] = ["0"]
        transcript.meta["contains"] = [transcript.id]
        transcript.meta["contains_count"] = 1
    
    def add_dists_to_ends(self, tm, contains):
        tm.meta["3p_dists_to_3p"].append(
            str(contains.TSS - tm.TSS)
        )
        
        tm.meta["5p_dists_to_5p"].append(
            str(abs(contains.TES - tm.TES))
        )

    def update_contains(self, tm, contains):
        tm.meta["contains"].extend(contains.meta["contains"])
        tm.meta["contains"].append(contains.id)
        tm.meta["contains_count"] += contains.meta["contains_count"]