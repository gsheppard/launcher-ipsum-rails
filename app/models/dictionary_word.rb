class DictionaryWord < ActiveRecord::Base
  def self.get_paragraph
    sent_per_para = rand(4..12)
    lines = DictionaryWord.last.id
    sent_ctr = 0
    paragraph = []

    while sent_ctr < sent_per_para
      words_per_sent = rand(5..10)
      word_ctr = 0
      sentence = []

      while word_ctr < words_per_sent
        sentence << DictionaryWord.find(rand(2...lines)).word
        word_ctr += 1
      end

      sentence.first.capitalize!
      paragraph << sentence
      sent_ctr += 1

    end

    paragraph
  end

  def self.build_paragraphs num_paragraphs
    ipsum = []
    ctr = 0
    while ctr < num_paragraphs
      ipsum << get_paragraph
      ctr += 1
    end
    ipsum
  end
end
