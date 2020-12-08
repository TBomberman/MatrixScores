class ScoreItem:
    def __init__(self,item_score,item_string):
        self.item_score= item_score
        self.item_string = item_string

    def __hash__(self):
        return hash((self.item_score, self.item_string))

    def __eq__(self, other):
        return (self.item_score, self.item_string) == (other.item_score, other.item_string)

    def __lt__(self, other):
        if self.item_score < other.item_score:
            return True
        return self.item_string < other.item_string