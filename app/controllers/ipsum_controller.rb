class IpsumController < ApplicationController

  def index

  end

  def post
    if input <= 0 || input > 20
      redirect_to root_path, alert: 'Please enter a valid number'
    else
      @generated_ipsum = DictionaryWord.build_paragraphs(input)
      render :index
    end
  end

  private
  def input
    params[:ipsum][:num].to_i
  end

end
