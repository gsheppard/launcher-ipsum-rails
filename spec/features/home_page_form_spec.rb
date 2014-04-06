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

  scenario 'fill in form with valid info' do
    visit root_path

    fill_in 'How many paragraphs?', with: '5'
    click_button 'Launch!'

    expect(page.all('p.ipsum').count).to eq(5)
  end

  scenario 'fill in form with blank info' do
    visit root_path

    click_button 'Launch!'

    expect(page.all('p.ipsum').count).to eq(0)
    expect(page).to have_content('Please enter a valid number')
  end

  scenario 'fill in with a number too large' do
    visit root_path

    fill_in 'How many paragraphs?', with: '21'
    click_button 'Launch!'

    expect(page.all('p.ipsum').count).to eq(0)
    expect(page).to have_content('Please enter a valid number')
  end

  scenario 'fill in with letters' do
    visit root_path

    fill_in 'How many paragraphs?', with: 'abc'
    click_button 'Launch!'

    expect(page.all('p.ipsum').count).to eq(0)
    expect(page).to have_content('Please enter a valid number')
  end


end
