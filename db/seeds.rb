# This file should contain all the record creation needed to seed the database with its default values.
# The data can then be loaded with the rake db:seed (or created alongside the db with db:setup).
#
# Examples:
#
#   cities = City.create([{ name: 'Chicago' }, { name: 'Copenhagen' }])
#   Mayor.create(name: 'Emanuel', city: cities.first)

progressbar = ProgressBar.create(total: nil)
File.open(File.join(Rails.root, 'db', 'dictionary.md')).each_line do |row|
  DictionaryWord.find_or_create_by(word: row.chomp)
  progressbar.increment
end
