require 'spec_helper'

feature 'home page with form', %q{
  As a user browsing
  I want to see a form on the homepage
  So I can generate that sweet ipsum
} do

  scenario 'see form on homepage' do
    visit root_path

    expect(page).to have_content('How many paragraphs?')
    expect(page).to have_button('Launch!')
  end

  scenario 'fill in form with valid info'

end
