Launchipsum::Application.routes.draw do

  root to: 'ipsum#index'
  post '/', to: 'ipsum#post'

end
