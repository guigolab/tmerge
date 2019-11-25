from models.transcript import Transcript

def build(exons):
    """
    Build a list of transcripts.

    Warning! Input exon data is mutated.
    """
    transcripts = []

    """
    Starting at exons[0]:
        - Build a list of all exons in exons list that share the same transcript_id as exon [0]
        - Remove each exon in generated list from exons list
        - Create a new Transcript from this list
        - Append Transcript to transcripts
        - Continue at next exon in exons list

    Note this may not be the fastest way to achieve this.
    TODO: find a faster method of building list of transcripts
    """
    while exons:
        from_same_transcript = [e for e in exons if exons[0].transcript_id == e.transcript_id]
        transcripts.append(Transcript(from_same_transcript))
        for e in from_same_transcript:
            exons.remove(e)
    
    return transcripts